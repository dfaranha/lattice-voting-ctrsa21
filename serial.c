#include <string.h>

#include "serial.h"
#include "commit.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/* Unary runs longer than this escape to a full-width value. */
#define RICE_ESCAPE 24

static int bits_for(uint64_t range) {
	int n = 0;

	while (range > 0) {
		n++;
		range >>= 1;
	}
	return n;
}

static void put_bits(bitwriter_t *w, uint64_t v, int n) {
	for (int i = 0; i < n; i++) {
		size_t at = w->bits + i;

		if (at / 8 >= w->cap) {
			w->overflow = 1;
			return;
		}
		if ((v >> i) & 1) {
			w->buf[at / 8] |= (uint8_t) (1u << (at % 8));
		} else {
			w->buf[at / 8] &= (uint8_t) ~(1u << (at % 8));
		}
	}
	w->bits += n;
}

static uint64_t get_bits(bitreader_t *r, int n) {
	uint64_t v = 0;

	for (int i = 0; i < n; i++) {
		size_t at = r->bits + i;

		if (at / 8 >= r->cap) {
			r->overflow = 1;
			return 0;
		}
		if ((r->buf[at / 8] >> (at % 8)) & 1) {
			v |= (uint64_t) 1 << i;
		}
	}
	r->bits += n;
	return v;
}

/* Zigzag, so that a small negative is a small unsigned. */
static uint64_t zigzag(int64_t x) {
	return ((uint64_t) x << 1) ^ (uint64_t) (x >> 63);
}

static int64_t unzigzag(uint64_t u) {
	return (int64_t) (u >> 1) ^ -(int64_t) (u & 1);
}

/*
 * Golomb-Rice parameter for a discrete Gaussian of width sigma. The mean of
 * the zigzagged coefficient is 2 sigma sqrt(2/pi), and Rice is near-optimal
 * when the split falls there, so k is the nearest power of two below it.
 */
static int rice_k(ulong sigma) {
	uint64_t mean = (uint64_t) (1.5957691 * (double) sigma);
	int k = 0;

	while ((mean >> (k + 1)) != 0) {
		k++;
	}
	return k;
}

/*
 * A coefficient, as quotient in unary and k low bits. The unary run is capped:
 * beyond the cap the value is written at full width, so that a coefficient far
 * outside the bound costs bits instead of costing correctness. There is no
 * width for it to overflow, which is why this no longer refuses anything.
 */
static void put_rice(bitwriter_t *w, int64_t x, int k) {
	uint64_t u = zigzag(x);
	uint64_t q = u >> k;

	if (q >= RICE_ESCAPE) {
		for (int i = 0; i < RICE_ESCAPE; i++) {
			put_bits(w, 1, 1);
		}
		put_bits(w, u, 64);
		return;
	}
	for (uint64_t i = 0; i < q; i++) {
		put_bits(w, 1, 1);
	}
	put_bits(w, 0, 1);
	if (k > 0) {
		put_bits(w, u & (((uint64_t) 1 << k) - 1), k);
	}
}

static int64_t get_rice(bitreader_t *r, int k) {
	int q = 0;
	uint64_t u;

	while (q < RICE_ESCAPE && get_bits(r, 1) == 1) {
		q++;
	}
	if (q >= RICE_ESCAPE) {
		return unzigzag(get_bits(r, 64));
	}
	u = (uint64_t) q << k;
	if (k > 0) {
		u |= get_bits(r, k);
	}
	return unzigzag(u);
}

/* Centre a residue into (-p/2, p/2]. */
static int64_t centre(ulong a) {
	return (a > MODP / 2) ? (int64_t) a - (int64_t) MODP : (int64_t) a;
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void serial_writer_init(bitwriter_t *w, uint8_t *buf, size_t cap) {
	w->buf = buf;
	w->cap = cap;
	w->bits = 0;
	w->overflow = 0;
	memset(buf, 0, cap);
}

void serial_reader_init(bitreader_t *r, const uint8_t *buf, size_t cap) {
	r->buf = buf;
	r->cap = cap;
	r->bits = 0;
	r->overflow = 0;
}

int serial_bits_uniform(void) {
	return bits_for(MODP - 1);
}

int serial_bits_gauss(ulong sigma) {
	/* The width a flat encoding would need: six standard deviations either
	 * side, so 12 sigma values in all. Nothing writes at this width any more
	 * -- serial_put_gauss codes the distribution instead -- but it is what the
	 * saving is measured against. */
	return bits_for(12 * sigma);
}

size_t serial_bytes(const bitwriter_t *w) {
	return (w->bits + 7) / 8;
}

void serial_put_uniform(bitwriter_t *w, nmod_poly_t a[2]) {
	int n = serial_bits_uniform();

	for (int k = 0; k < NCRT; k++) {
		for (int i = 0; i < DEGCRT; i++) {
			put_bits(w, nmod_poly_get_coeff_ui(a[k], i), n);
		}
	}
}

void serial_get_uniform(bitreader_t *r, nmod_poly_t a[2]) {
	int n = serial_bits_uniform();

	for (int k = 0; k < NCRT; k++) {
		nmod_poly_zero(a[k]);
		nmod_poly_fit_length(a[k], DEGCRT);
		for (int i = 0; i < DEGCRT; i++) {
			nmod_poly_set_coeff_ui(a[k], i, get_bits(r, n));
		}
	}
}

void serial_put_gauss(bitwriter_t *w, nmod_poly_t a[2], ulong sigma) {
	int k = rice_k(sigma);
	nmod_poly_t rec;

	nmod_poly_init(rec, MODP);
	pcrt_poly_rec(rec, a);
	for (int i = 0; i < DEGREE; i++) {
		put_rice(w, centre(nmod_poly_get_coeff_ui(rec, i)), k);
	}
	nmod_poly_clear(rec);
}

void serial_get_gauss(bitreader_t *r, nmod_poly_t a[2], ulong sigma) {
	int k = rice_k(sigma);
	nmod_poly_t rec;

	nmod_poly_init(rec, MODP);
	nmod_poly_fit_length(rec, DEGREE);
	for (int i = 0; i < DEGREE; i++) {
		int64_t c = get_rice(r, k);
		ulong v = (c < 0) ? MODP - ((ulong) (-c) % MODP) : (ulong) c % MODP;

		nmod_poly_set_coeff_ui(rec, i, v == MODP ? 0 : v);
	}
	pcrt_poly_reduce(a[0], rec, 0);
	pcrt_poly_reduce(a[1], rec, 1);
	nmod_poly_clear(rec);
}


void serial_put_byte(bitwriter_t *w, uint8_t v) {
	put_bits(w, v, 8);
}

uint8_t serial_get_byte(bitreader_t *r) {
	return (uint8_t) get_bits(r, 8);
}
