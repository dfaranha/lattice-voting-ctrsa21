#include <string.h>

#include "serial.h"
#include "commit.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

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
	/* Six standard deviations either side, so 12 sigma values in all. */
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
	int n = serial_bits_gauss(sigma);
	int64_t half = 6 * (int64_t) sigma;
	nmod_poly_t rec;

	nmod_poly_init(rec, MODP);
	pcrt_poly_rec(rec, a);
	for (int i = 0; i < DEGREE; i++) {
		int64_t c = centre(nmod_poly_get_coeff_ui(rec, i));

		/* The modelling assumption is that six standard deviations is enough.
		 * Refuse rather than truncate, so that a run needing more says so
		 * instead of quietly reporting a size that could not be decoded. The
		 * flag rather than an assertion, so that it holds under NDEBUG too. */
		if (c < -half || c > half) {
			w->overflow = 1;
			nmod_poly_clear(rec);
			return;
		}
		put_bits(w, (uint64_t) (c + half), n);
	}
	nmod_poly_clear(rec);
}

void serial_get_gauss(bitreader_t *r, nmod_poly_t a[2], ulong sigma) {
	int n = serial_bits_gauss(sigma);
	int64_t half = 6 * (int64_t) sigma;
	nmod_poly_t rec;

	nmod_poly_init(rec, MODP);
	nmod_poly_fit_length(rec, DEGREE);
	for (int i = 0; i < DEGREE; i++) {
		int64_t c = (int64_t) get_bits(r, n) - half;
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
