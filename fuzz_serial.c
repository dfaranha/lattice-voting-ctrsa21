/*
 * libFuzzer harness for the proof decoder in serial.c.
 *
 * The decoder reads what an untrusted prover sends, so it must not misbehave on
 * any input, which the sanitizers check, and it must not accept two byte
 * strings for one proof, which is the oracle here: an input the decoder
 * accepts must re-encode to exactly the same bytes.
 *
 * The first byte is a plan: its top two bits give the number of elements, one
 * to four, and each pair of the rest the kind of one element, with the fourth
 * element's kind taken from the second byte. The rest is the encoding.
 *
 *   make fuzz_serial
 *   ./fuzz_serial -seeds fuzz-corpus     # write a seed corpus of valid encodings
 *   ./fuzz_serial -max_len=65536 fuzz-corpus
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "serial.h"
#include "commit.h"
#include "lnp.h"

#define MAX_ELEMS 4

enum { KIND_UNIFORM, KIND_GAUSS, KIND_PROJ, KIND_BYTE };

/* On this branch the third kind is the published projection, PROJ centred
 * residues at SIGMA_P, which has a codec of its own. */

static int plan_parse(const uint8_t *data, size_t size, int kinds[MAX_ELEMS],
		size_t *used) {
	int n;

	if (size < 2) {
		return 0;
	}
	n = 1 + (data[0] >> 6);
	for (int i = 0; i < 3; i++) {
		kinds[i] = (data[0] >> (2 * i)) & 3;
	}
	kinds[3] = data[1] & 3;
	*used = 2;
	return n;
}

static void walk(bitwriter_t *w, bitreader_t *r, int kind, nmod_poly_t a[2],
		uint8_t *byte, ulong proj[PROJ]) {
	switch (kind) {
		case KIND_UNIFORM:
			if (w != NULL) {
				serial_put_uniform(w, a);
			} else {
				serial_get_uniform(r, a);
			}
			break;
		case KIND_GAUSS:
			if (w != NULL) {
				serial_put_gauss(w, a, SIGMA_C);
			} else {
				serial_get_gauss(r, a, SIGMA_C);
			}
			break;
		case KIND_PROJ:
			if (w != NULL) {
				serial_put_proj(w, proj, PROJ, SIGMA_P);
			} else {
				serial_get_proj(r, proj, PROJ, SIGMA_P);
			}
			break;
		default:
			if (w != NULL) {
				serial_put_byte(w, *byte);
			} else {
				*byte = serial_get_byte(r);
			}
			break;
	}
}

static int setup(void) {
	static int done = 0;

	if (!done) {
		commit_setup();
		done = 1;
	}
	return 0;
}

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
	int kinds[MAX_ELEMS], n;
	size_t used;
	nmod_poly_t a[MAX_ELEMS][2];
	uint8_t bytes[MAX_ELEMS] = { 0 };
	static ulong proj[MAX_ELEMS][PROJ];
	bitreader_t r;

	setup();
	n = plan_parse(data, size, kinds, &used);
	if (n == 0) {
		return 0;
	}
	for (int i = 0; i < n; i++) {
		nmod_poly_init(a[i][0], MODP);
		nmod_poly_init(a[i][1], MODP);
	}

	serial_reader_init(&r, data + used, size - used);
	for (int i = 0; i < n; i++) {
		walk(NULL, &r, kinds[i], a[i], &bytes[i], proj[i]);
	}
	if (serial_reader_finish(&r)) {
		size_t len = size - used;
		uint8_t *buf = (uint8_t *) malloc(len > 0 ? len : 1);
		bitwriter_t w;

		serial_writer_init(&w, buf, len);
		for (int i = 0; i < n; i++) {
			walk(&w, NULL, kinds[i], a[i], &bytes[i], proj[i]);
		}
		if (w.overflow || serial_bytes(&w) != len
				|| memcmp(buf, data + used, len) != 0) {
			fprintf(stderr, "accepted an encoding that re-encodes differently\n");
			abort();
		}
		free(buf);
	}

	for (int i = 0; i < n; i++) {
		nmod_poly_clear(a[i][0]);
		nmod_poly_clear(a[i][1]);
	}
	return 0;
}

/* Writes valid encodings of every plan, so that the fuzzer starts from inputs
 * the decoder accepts rather than having to find one. Called from
 * LLVMFuzzerInitialize when the first argument is -seeds. */
static void write_seeds(const char *dir) {
	static uint8_t buf[1 << 18];
	flint_rand_t rand;
	int count = 0;

	setup();
	flint_rand_init(rand);
	for (int plan = 0; plan < 256; plan += 5) {
		for (int variant = 0; variant < 3; variant++) {
			int kinds[MAX_ELEMS], n;
			size_t used;
			uint8_t head[2] = { (uint8_t) plan, (uint8_t) (plan * 7 + variant) };
			nmod_poly_t a[2], t;
			uint8_t byte = (uint8_t) (plan ^ variant);
			static ulong proj[PROJ];
			bitwriter_t w;
			char path[512];
			FILE *f;

			n = plan_parse(head, 2, kinds, &used);
			nmod_poly_init(a[0], MODP);
			nmod_poly_init(a[1], MODP);
			nmod_poly_init(t, MODP);
			serial_writer_init(&w, buf, sizeof(buf));
			for (int i = 0; i < n; i++) {
				nmod_poly_zero(t);
				for (int j = 0; j < DEGREE; j++) {
					/* Zero, uniform, or small with an occasional outlier that
					 * forces an escape. */
					ulong v = 0;

					if (variant == 1 || kinds[i] == KIND_UNIFORM) {
						v = n_randint(rand, MODP);
					} else if (variant == 2) {
						long c = (long) n_randint(rand, 4 * SIGMA_C) - 2 * SIGMA_C;

						if (n_randint(rand, 64) == 0) {
							c *= 40;
						}
						v = (c < 0) ? (ulong) (MODP + c) : (ulong) c;
					}
					nmod_poly_set_coeff_ui(t, j, v);
				}
				pcrt_poly_reduce(a[0], t, 0);
				pcrt_poly_reduce(a[1], t, 1);
				for (int j = 0; j < PROJ; j++) {
					long c = (variant == 0) ? 0 :
							(long) n_randint(rand, 8 * SIGMA_P) - 4 * SIGMA_P;

					proj[j] = (c < 0) ? (ulong) (MODP + c) : (ulong) c;
				}
				walk(&w, NULL, kinds[i], a, &byte, proj);
			}
			nmod_poly_clear(a[0]);
			nmod_poly_clear(a[1]);
			nmod_poly_clear(t);
			if (w.overflow) {
				continue;
			}
			snprintf(path, sizeof(path), "%s/seed-%03d-%d", dir, plan, variant);
			f = fopen(path, "wb");
			if (f == NULL) {
				perror(path);
				exit(1);
			}
			fwrite(head, 1, 2, f);
			fwrite(buf, 1, serial_bytes(&w), f);
			fclose(f);
			count++;
		}
	}
	flint_rand_clear(rand);
	printf("wrote %d seeds to %s\n", count, dir);
}

extern "C" int LLVMFuzzerInitialize(int *argc, char ***argv) {
	if (*argc >= 3 && strcmp((*argv)[1], "-seeds") == 0) {
		write_seeds((*argv)[2]);
		exit(0);
	}
	return 0;
}
