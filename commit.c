/**
 * @file
 *
 * Implementation of the lattice-based commitment scheme.
 *
 * @ingroup commit
 */

#include <math.h>
#include <flint/nmod.h>

#include "param.h"
#include "commit.h"
#include "test.h"
#include "bench.h"
#include "gaussian.h"
#include "randombytes.h"
#include "fastrandombytes.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/* The first square root of -1. */
#define	P0		499810066026
/* The second square root of -1. */
#define P1		599701561891

/* Polynomial defining the cyclotomic ring. */
static nmod_poly_t cyclo_poly;

/* Pair of irreducible polynomials for CRT representation. */
static pcrt_poly_t irred;

/* Inverses of the irreducible polynomials for CRT reconstruction. */
static pcrt_poly_t inv;

/* Scratch space for the multiplication routines. Like the rest of this module,
 * they are not reentrant. */
static nmod_poly_t mul_tmp;

/**
 * Reduce a polynomial modulo (x^m + r).
 *
 * Both moduli used here have this shape: the cyclotomic polynomial is
 * x^DEGREE + 1, and each CRT factor is x^DEGCRT + r_i. Reduction is then just a
 * fold of each successive block of m coefficients onto the lowest one, scaled
 * by (-r)^k, which is far cheaper than the generic polynomial division that
 * nmod_poly_mulmod performs. Measured on the parameters of this code, it makes
 * a ring multiplication about 3x faster.
 *
 * The input may have any degree. The output may alias the input: the copy loop
 * is then a no-op and the fold loop only reads coefficients at or above m,
 * which it never writes.
 *
 * @param[out] c		- the reduced polynomial.
 * @param[in] a			- the polynomial to reduce.
 * @param[in] m			- the degree of the modulus.
 * @param[in] r			- the constant coefficient of the modulus.
 */
static void fold_mod(nmod_poly_t c, const nmod_poly_t a, slong m, ulong r) {
	nmod_t mod = a->mod;
	slong len = a->length;
	ulong f = 1;

	nmod_poly_fit_length(c, m);
	for (slong j = 0; j < m; j++) {
		c->coeffs[j] = (j < len) ? a->coeffs[j] : 0;
	}
	/* The block of coefficients at x^(k*m) contributes a factor (-r)^k. */
	for (slong base = m; base < len; base += m) {
		f = nmod_neg(nmod_mul(f, r, mod), mod);
		for (slong j = 0; j < m && base + j < len; j++) {
			c->coeffs[j] = nmod_add(c->coeffs[j],
					nmod_mul(a->coeffs[base + j], f, mod), mod);
		}
	}
	c->length = m;
	_nmod_poly_normalise(c);
}

/**
 * Test if the l2-norm is within bounds (4 * sigma * sqrt(N)).
 *
 * @param[in] r 			- the polynomial to compute the l2-norm.
 * @return the computed norm.
 */
static int test_norm(nmod_poly_t r) {
	/* The width the masks are sampled at, and so the one an extracted opening
	 * is bounded by. The binding estimates use the same SIGMA_C, through the
	 * extractable bound 16 * SIGMA_C * sqrt(nu * N). The paper's lower bound
	 * 11 * nu * beta * sqrt(k * N) on that width, which this used to take,
	 * accepted only openings 2.1 times shorter than extraction yields. */
	const uint64_t sigma_sqr = (uint64_t) SIGMA_C * SIGMA_C;

	/* Compare the squared norm to (4 * sigma * sqrt(N))^2 = 16 * sigma^2 * N.
	 * The opening is the prover's input, so this has to be the check that
	 * cannot overflow: summing the squares in a word wraps for a handful of
	 * coefficients near MODP / 2, which let an enormous opening through. */
	return commit_norm2_leq(r, (uint64_t) 16 * sigma_sqr * DEGREE);
}

/**
 * Test if a polynomial is a valid opening challenge: a non-zero difference
 * of two challenges, so that ||f||_oo <= 2 and ||f||_1 <= 2 * NONZERO.
 *
 * An opening (m, r, f) proves f * c = B * r + f * m, which says nothing about m
 * when f is zero: r = 0 and f = 0 open any commitment to any message. A short
 * non-zero f is invertible by Lemma 1, which is what makes the opening bind.
 *
 * @param[in] f 			- the challenge in coefficient representation.
 * @return 1 if it is a valid opening challenge, 0 otherwise.
 */
static int test_chall_diff(nmod_poly_t f) {
	int64_t coeff;
	uint64_t l1 = 0;

	if (nmod_poly_is_zero(f)) {
		return 0;
	}
	for (int i = 0; i < DEGREE; i++) {
		coeff = nmod_poly_get_coeff_ui(f, i);
		if (coeff > MODP / 2) {
			coeff -= MODP;
		}
		if (coeff < -2 || coeff > 2) {
			return 0;
		}
		l1 += (coeff < 0 ? -coeff : coeff);
	}
	return l1 <= 2 * NONZERO;
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

/**
 * Sample a uniformly random double in [0, 1) with 53 bits of precision.
 *
 * Rejection sampling must not reuse the Fiat-Shamir stream, so this draws from
 * the operating system rather than from fastrandombytes.
 *
 * @return the sampled value.
 */
double commit_uniform_double(void) {
	uint64_t bits;

	randombytes((unsigned char *)&bits, sizeof(bits));
	/* Keep the top 53 bits, the precision of the mantissa of a double. */
	return (double)(bits >> 11) * 0x1.0p-53;
}

/* z and v are in coefficient representation: the test needs their centred
 * coefficients, and the provers keep them that way so as not to reconstruct
 * them from CRT form on every attempt. */
void commit_rej_accumulate(int64_t *dot, int64_t *norm, nmod_poly_t z[],
		nmod_poly_t v[], int width) {
	int64_t c0, c1;

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < DEGREE; j++) {
			c0 = nmod_poly_get_coeff_ui(z[i], j);
			c1 = nmod_poly_get_coeff_ui(v[i], j);
			if (c0 > MODP / 2)
				c0 -= MODP;
			if (c1 > MODP / 2)
				c1 -= MODP;
			*dot += c0 * c1;
			*norm += c1 * c1;
		}
	}
}

void commit_rej_respond(nmod_poly_t z[][2], nmod_poly_t y[], nmod_poly_t v[],
		nmod_poly_t d, int width, int64_t *dot, int64_t *norm) {
	nmod_poly_t dv;

	nmod_poly_init(dv, MODP);
	for (int i = 0; i < width; i++) {
		commit_poly_mul_sparse(dv, d, v[i]);
		nmod_poly_add(y[i], y[i], dv);
		commit_rej_accumulate(dot, norm, &y[i], &dv, 1);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_reduce(z[i][k], y[i], k);
		}
	}
	nmod_poly_clear(dv);
}

int commit_rej_decide(int64_t dot, int64_t norm, uint64_t s2) {
	double r = (-2.0 * (double) dot + (double) norm) / (2.0 * (double) s2);

	return commit_uniform_double() > exp(r) / 1.75;
}

int commit_rej_sampling(nmod_poly_t z[], nmod_poly_t v[], uint64_t s2,
		int width) {
	int64_t dot = 0, norm = 0;

	commit_rej_accumulate(&dot, &norm, z, v, width);
	return commit_rej_decide(dot, norm, s2);
}

// Initialize commitment scheme.
void commit_setup() {
	nmod_poly_init(cyclo_poly, MODP);
	nmod_poly_init(mul_tmp, MODP);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(irred[i], MODP);
		nmod_poly_init(inv[i], MODP);
	}

	// Initialize polynomial as x^N + 1. */
	nmod_poly_set_coeff_ui(cyclo_poly, DEGREE, 1);
	nmod_poly_set_coeff_ui(cyclo_poly, 0, 1);

	// Initialize two factors of the polynomial for CRT representation.
	nmod_poly_set_coeff_ui(irred[0], DEGCRT, 1);
	nmod_poly_set_coeff_ui(irred[0], 0, P0);
	nmod_poly_set_coeff_ui(irred[1], DEGCRT, 1);
	nmod_poly_set_coeff_ui(irred[1], 0, P1);

	nmod_poly_invmod(inv[0], irred[0], irred[1]);
	nmod_poly_invmod(inv[1], irred[1], irred[0]);
	nmod_poly_mul(inv[1], inv[1], irred[1]);
}

// Finalize commitment scheme.
void commit_finish() {
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(irred[i]);
		nmod_poly_clear(inv[i]);
	}
	nmod_poly_clear(cyclo_poly);
	nmod_poly_clear(mul_tmp);
}

// Return polynomial defining Rp.
nmod_poly_t *commit_poly() {
	return &cyclo_poly;
}

// Return irreducible polynomials defining CRT representation.
nmod_poly_t *commit_irred(int i) {
	return &irred[i];
}

// Multiply two polynomials modulo the i-th CRT factor.
void pcrt_poly_mulmod(nmod_poly_t c, const nmod_poly_t a, const nmod_poly_t b,
		int i) {
	nmod_poly_mul(mul_tmp, a, b);
	fold_mod(c, mul_tmp, DEGCRT, nmod_poly_get_coeff_ui(irred[i], 0));
}

// Multiply two polynomials in Rp.
void commit_poly_mulmod(nmod_poly_t c, const nmod_poly_t a,
		const nmod_poly_t b) {
	nmod_poly_mul(mul_tmp, a, b);
	fold_mod(c, mul_tmp, DEGREE, 1);
}

// Multiply by a sparse polynomial in Rp.
void commit_poly_mul_sparse(nmod_poly_t c, const nmod_poly_t d,
		const nmod_poly_t a) {
	nmod_t mod = a->mod;
	slong la = FLINT_MIN(a->length, DEGREE);
	ulong *out = (ulong *) flint_calloc(DEGREE, sizeof(ulong));

	/* Each non-zero coefficient of d adds a signed rotation of a: x^p * x^j
	 * is x^(p + j), or -x^(p + j - DEGREE) once it wraps past x^DEGREE = -1.
	 * The challenges have NONZERO coefficients of +-1, so this costs
	 * NONZERO * DEGREE additions where a full product costs a multiplication
	 * at full degree. */
	for (slong p = 0; p < FLINT_MIN(d->length, DEGREE); p++) {
		ulong v = d->coeffs[p];

		if (v == 0) {
			continue;
		}
		for (slong j = 0; j < la; j++) {
			ulong t = a->coeffs[j];

			if (v != 1) {
				t = (v == mod.n - 1) ? nmod_neg(t, mod) : nmod_mul(t, v, mod);
			}
			if (p + j < DEGREE) {
				out[p + j] = nmod_add(out[p + j], t, mod);
			} else {
				out[p + j - DEGREE] = nmod_sub(out[p + j - DEGREE], t, mod);
			}
		}
	}
	nmod_poly_fit_length(c, DEGREE);
	for (slong j = 0; j < DEGREE; j++) {
		c->coeffs[j] = out[j];
	}
	c->length = DEGREE;
	_nmod_poly_normalise(c);
	flint_free(out);
}

// Reduce a polynomial into the i-th CRT component.
void pcrt_poly_reduce(nmod_poly_t c, const nmod_poly_t a, int i) {
	fold_mod(c, a, DEGCRT, nmod_poly_get_coeff_ui(irred[i], 0));
}

// Recover polynomial from CRT representation.
void pcrt_poly_rec(nmod_poly_t c, pcrt_poly_t a) {
	nmod_poly_sub(c, a[0], a[1]);
	nmod_poly_mul(mul_tmp, c, inv[1]);
	/* a[1] is already reduced, so it can be added after the fold. */
	fold_mod(c, mul_tmp, DEGREE, 1);
	nmod_poly_add(c, c, a[1]);
}

// Compute squared l2-norm.
/* Only meaningful for short polynomials: the accumulator overflows for
 * coefficients near MODP. Use commit_norm2_leq to check a bound on input that
 * a malicious party could have chosen. */
uint64_t commit_norm2_sqr(nmod_poly_t r) {
	int64_t coeff;
	uint64_t norm = 0;

	/* Compute norm^2. Unsigned, so that an overflow wraps rather than being
	 * undefined, which it was for the signed accumulator this used to be. */
	for (int i = 0; i < DEGREE; i++) {
		coeff = nmod_poly_get_coeff_ui(r, i);
		if (coeff > MODP / 2)
			coeff -= MODP;
		norm += (uint64_t) (coeff * coeff);
	}
	return norm;
}

// Compute the l-infinity norm.
uint64_t commit_norm_inf(nmod_poly_t r) {
	int64_t coeff;
	uint64_t max = 0;

	for (int i = 0; i < DEGREE; i++) {
		coeff = nmod_poly_get_coeff_ui(r, i);
		if (coeff > MODP / 2) {
			coeff -= MODP;
		}
		if (coeff < 0) {
			coeff = -coeff;
		}
		if ((uint64_t) coeff > max) {
			max = coeff;
		}
	}
	return max;
}

// Test whether the squared l2-norm is at most a bound, without overflowing.
int commit_norm2_leq(nmod_poly_t r, uint64_t bound) {
	int64_t coeff;
	uint64_t norm = 0;

	for (int i = 0; i < DEGREE; i++) {
		coeff = nmod_poly_get_coeff_ui(r, i);
		if (coeff > MODP / 2) {
			coeff -= MODP;
		}
		if (coeff < 0) {
			coeff = -coeff;
		}
		/* Bail out before squaring could overflow the accumulator. Comparing
		 * against bound alone is not enough: bound is a squared norm, so a
		 * coefficient below it can still square past 2^64. Dividing instead
		 * of multiplying decides the same question without overflowing. This
		 * was unreachable while MODP stayed under 2^32, since no centred
		 * coefficient could then reach the square root of the accumulator. */
		if (coeff != 0 && (uint64_t) coeff > bound / (uint64_t) coeff) {
			return 0;
		}
		norm += (uint64_t) coeff * coeff;
		if (norm > bound) {
			return 0;
		}
	}
	return 1;
}

// Initialise a commitment key pair.
void commit_keyinit(commitkey_t *key) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(key->B1[i][j][k], MODP);
			}
		}
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_init(key->b2[i][j], MODP);
		}
	}
}

/* Coefficients are absorbed at the narrowest width that holds MODP, so that
 * the same code serves a modulus below 2^32, where four bytes suffice, and one
 * above it, such as the 40-bit prime here, which needs eight. */
#if MODP < (1ULL << 32)
typedef uint32_t hash_coeff_t;
#else
typedef uint64_t hash_coeff_t;
#endif

static_assert(MODP - 1 <= (hash_coeff_t) ~(hash_coeff_t) 0,
		"commit_hash_poly absorbs coefficients too narrow to hold MODP");

static void hash_coeffs(SHA256Context *sha, nmod_poly_t p, int len) {
	hash_coeff_t buf[DEGREE];

	for (int i = 0; i < len; i++) {
		buf[i] = (hash_coeff_t) nmod_poly_get_coeff_ui(p, i);
	}
	SHA256Input(sha, (const uint8_t *)buf, len * sizeof(hash_coeff_t));
}

void commit_hash_poly(SHA256Context *sha, nmod_poly_t p) {
	hash_coeffs(sha, p, DEGREE);
}

void commit_hash_crt(SHA256Context *sha, nmod_poly_t p) {
	if (nmod_poly_length(p) > DEGCRT) {
		fprintf(stderr, "commit_hash_crt: CRT component is not reduced\n");
		abort();
	}
	hash_coeffs(sha, p, DEGCRT);
}

// Generate a key pair.
void commit_keygen(commitkey_t *key) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_zero(key->B1[i][j][k]);
				if (i == j) {
					nmod_poly_set_coeff_ui(key->B1[i][j][k], 0, 1);
				}
			}
		}
	}

	for (int i = 0; i < HEIGHT; i++) {
		for (int j = HEIGHT; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				commit_sample_rand(key->B1[i][j][k], NULL, DEGCRT);
			}
		}
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_zero(key->b2[i][j]);
			if (i == HEIGHT) {
				nmod_poly_set_coeff_ui(key->b2[i][j], 0, 1);
			}
			if (i > HEIGHT) {
				commit_sample_rand(key->b2[i][j], NULL, DEGCRT);
			}
		}
	}

	SHA256Context sha;

	SHA256Reset(&sha);
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				commit_hash_crt(&sha, key->B1[i][j][k]);
			}
		}
	}
	for (int j = 0; j < WIDTH; j++) {
		for (int k = 0; k < NCRT; k++) {
			commit_hash_crt(&sha, key->b2[j][k]);
		}
	}
	SHA256Result(&sha, key->digest);
}

// Free a commitment key pair.
void commit_keyfree(commitkey_t *key) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(key->B1[i][j][k]);
			}
		}
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_clear(key->b2[i][j]);
		}
	}
}

// Sample a short polynomial.
void commit_sample_short(nmod_poly_t r) {
	uint64_t buf;
	/* This holds MODP - 1 + d for d in {0, 1, 2}, so it has to be as wide as
	 * the modulus. It was uint32_t, which silently truncated as soon as MODP
	 * grew past 2^32 and turned the ternary randomness into coefficients near
	 * MODP mod 2^32, at which point rejection sampling never accepts. */
	ulong coeff;
	int i, j, s;

	nmod_poly_zero(r);
	nmod_poly_fit_length(r, DEGREE);
	i = 0;
	s = 8 * sizeof(buf);
	j = s;

	do {
		if (j == s) {
			randombytes((unsigned char *)&buf, sizeof(buf));
			j = 0;
		}

		if (((buf >> j) & 1) & ((buf >> (j + 1)) & 1)) {
			j += 2;
		} else {
			coeff = MODP - 1 + ((buf >> j) & 3);
			nmod_poly_set_coeff_ui(r, i, coeff);
			i++;
			j += 2;
		}
	} while (i < DEGREE);
}

// Sample a short polynomial in CRT representation.
void commit_sample_short_crt(pcrt_poly_t r) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	commit_sample_short(t);
	for (int j = 0; j < NCRT; j++) {
		pcrt_poly_reduce(r[j], t, j);
	}
	nmod_poly_clear(t);
}

// Sample a uniformly random polynomial.
void commit_sample_rand(nmod_poly_t r, prg_t *prg, int degree) {
	/* 2^64 mod MODP. Words at or above 2^64 minus this fall in an incomplete
	 * copy of [0, MODP), so rejecting them leaves the reduction unbiased. */
	const uint64_t rem = (UINT64_MAX % MODP + 1) % MODP;
	uint64_t buf;

	nmod_poly_fit_length(r, degree);
	for (int i = 0; i < degree; i++) {
		do {
			if (prg == NULL) {
				fastrandombytes((unsigned char *)&buf, sizeof(buf));
			} else {
				prg_bytes(prg, (unsigned char *)&buf, sizeof(buf));
			}
		} while (rem != 0 && buf >= (uint64_t) 0 - rem);
		r->coeffs[i] = buf % MODP;
	}
	r->length = degree;
	_nmod_poly_normalise(r);
}

// Sample a random polynomial in CRT representation.
void commit_sample_rand_crt(pcrt_poly_t r, prg_t *prg) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	commit_sample_rand(t, prg, DEGREE);
	for (int i = 0; i < NCRT; i++) {
		pcrt_poly_reduce(r[i], t, i);
	}
	nmod_poly_clear(t);
}

// Sample a challenge.
void commit_sample_chall(nmod_poly_t f) {
	nmod_poly_zero(f);
	nmod_poly_t c[2];
	uint32_t buf;

	/* The two challenges of the difference c[0] - c[1], each with NONZERO
	 * coefficients of +-1. This 2 belongs to the challenge construction, not
	 * to the CRT decomposition. */
	for (int i = 0; i < 2; i++) {
		nmod_poly_init(c[i], MODP);
		nmod_poly_fit_length(c[i], DEGREE);
		for (int j = 0; j < NONZERO; j++) {
			uint32_t pos;

			do {
				randombytes((unsigned char *)&buf, sizeof(buf));
				pos = buf % DEGREE;
			} while (nmod_poly_get_coeff_ui(c[i], pos) != 0);
			nmod_poly_set_coeff_ui(c[i], pos, (buf >> 31) ? MODP - 1 : 1);
		}
	}
	nmod_poly_sub(f, c[0], c[1]);

	nmod_poly_clear(c[0]);
	nmod_poly_clear(c[1]);
}

// Sample a challenge in CRT representation.
void commit_sample_chall_crt(pcrt_poly_t f) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	commit_sample_chall(t);
	for (int i = 0; i < NCRT; i++) {
		pcrt_poly_reduce(f[i], t, i);
	}
	nmod_poly_clear(t);
}

// Sample a polynomial according to a Gaussian distribution.
void commit_sample_gauss(nmod_poly_t r) {
	int64_t coeff;
	for (int i = 0; i < DEGREE; i++) {
		coeff = discrete_gaussian(0.0);
		if (coeff < 0)
			coeff += MODP;
		nmod_poly_set_coeff_ui(r, i, coeff);
	}
}

// Sample a polynomial according to a Gaussian distribution in CRT rep.
void commit_sample_gauss_crt(nmod_poly_t r[2]) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	commit_sample_gauss(t);
	for (int i = 0; i < NCRT; i++) {
		pcrt_poly_reduce(r[i], t, i);
	}

	nmod_poly_clear(t);
}

void commit_sample_gauss_batch(nmod_poly_t t) {
	int64_t coeff;

	nmod_poly_zero(t);
	nmod_poly_fit_length(t, DEGREE);
	for (int i = 0; i < DEGREE; i++) {
		coeff = discrete_gaussian_batch(0.0);
		nmod_poly_set_coeff_ui(t, i,
				(coeff < 0) ? MODP - ((-coeff) % MODP) : (ulong)(coeff % MODP));
	}
}

void commit_sample_gauss_batch_crt(nmod_poly_t r[2]) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	commit_sample_gauss_batch(t);
	for (int i = 0; i < NCRT; i++) {
		pcrt_poly_reduce(r[i], t, i);
	}
	nmod_poly_clear(t);
}

// Initialise a commitment.
void commit_init(commit_t *com) {
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(com->c1[i], MODP);
		nmod_poly_init(com->c2[i], MODP);
	}
}

// Sample a polynomial according to a narrow Gaussian distribution in CRT rep.
void commit_sample_gauss_small_crt(nmod_poly_t r[NCRT]) {
	nmod_poly_t t;
	int64_t coeff;

	nmod_poly_init(t, MODP);
	for (int i = 0; i < DEGREE; i++) {
		coeff = discrete_gaussian_small(0.0);
		if (coeff < 0)
			coeff += MODP;
		nmod_poly_set_coeff_ui(t, i, coeff);
	}
	for (int i = 0; i < NCRT; i++) {
		pcrt_poly_reduce(r[i], t, i);
	}

	nmod_poly_clear(t);
}

// Commit to a message.
void commit_doit(commit_t *com, nmod_poly_t m, commitkey_t *key,
		pcrt_poly_t r[WIDTH]) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_zero(com->c1[i]);
		nmod_poly_zero(com->c2[i]);
	}

	// Compute B = [ B1 b2 ]^t * r_m.
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_mulmod(t, key->B1[i][j][k], r[j][k], k);
				nmod_poly_add(com->c1[k], com->c1[k], t);
				if (i == 0) {
					pcrt_poly_mulmod(t, key->b2[j][k], r[j][k], k);
					nmod_poly_add(com->c2[k], com->c2[k], t);
				}
			}
		}
	}

	// Convert m to CRT representation and accumulate.
	for (int i = 0; i < NCRT; i++) {
		pcrt_poly_reduce(t, m, i);
		nmod_poly_add(com->c2[i], com->c2[i], t);
	}

	nmod_poly_clear(t);
}

// Open a commitment on a message, randomness, factor.
int commit_open(commit_t *com, nmod_poly_t m, commitkey_t *key,
		pcrt_poly_t r[WIDTH], pcrt_poly_t f) {
	nmod_poly_t t;
	pcrt_poly_t c1, c2, _c1, _c2;
	int result = 0;

	nmod_poly_init(t, MODP);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(c1[i], MODP);
		nmod_poly_init(c2[i], MODP);
		nmod_poly_init(_c1[i], MODP);
		nmod_poly_init(_c2[i], MODP);
		nmod_poly_zero(c1[i]);
		nmod_poly_zero(c2[i]);
	}

	// Compute B = [ B1 b2 ]^t * r_m.
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_mulmod(t, key->B1[i][j][k], r[j][k], k);
				nmod_poly_add(c1[k], c1[k], t);
				if (i == 0) {
					pcrt_poly_mulmod(t, key->b2[j][k], r[j][k], k);
					nmod_poly_add(c2[k], c2[k], t);
				}
			}
		}
	}

	// Convert m to CRT representation before multiplication.
	for (int i = 0; i < NCRT; i++) {
		pcrt_poly_reduce(t, m, i);
		pcrt_poly_mulmod(t, t, f[i], i);
		nmod_poly_add(c2[i], c2[i], t);
	}

	for (int i = 0; i < NCRT; i++) {
		pcrt_poly_mulmod(_c1[i], com->c1[i], f[i], i);
		pcrt_poly_mulmod(_c2[i], com->c2[i], f[i], i);
	}

	/* Every component of r is bounded, not only the first three, which is
	 * all this used to check whatever the WIDTH. */
	pcrt_poly_rec(t, f);
	result = test_chall_diff(t);
	for (int i = 0; i < WIDTH; i++) {
		pcrt_poly_rec(t, r[i]);
		result &= test_norm(t);
	}
	for (int i = 0; i < NCRT; i++) {
		result &= nmod_poly_equal(_c1[i], c1[i]);
		result &= nmod_poly_equal(_c2[i], c2[i]);
	}

	nmod_poly_clear(t);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(c1[i]);
		nmod_poly_clear(c2[i]);
		nmod_poly_clear(_c1[i]);
		nmod_poly_clear(_c2[i]);
	}
	return result;
}

// Free a commitment.
void commit_free(commit_t *com) {
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(com->c1[i]);
		nmod_poly_clear(com->c2[i]);
	}
}

#ifdef MAIN
// Tests and benchmarks below.
static void test(flint_rand_t rand) {
	commitkey_t key;
	commit_t com, _com;
	nmod_poly_t m, _m, rho;
	pcrt_poly_t r[WIDTH], s[WIDTH], _s[WIDTH], f, _f;

	nmod_poly_init(m, MODP);
	nmod_poly_init(_m, MODP);
	nmod_poly_init(rho, MODP);
	for (int j = 0; j < NCRT; j++) {
		nmod_poly_init(f[j], MODP);
		nmod_poly_init(_f[j], MODP);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_init(r[i][j], MODP);
			nmod_poly_init(s[i][j], MODP);
			nmod_poly_init(_s[i][j], MODP);
		}
	}

	/* Generate a random message. */
	nmod_poly_randtest(m, rand, DEGREE);

	/* Generate commitment key. */
	commit_keyinit(&key);
	commit_keygen(&key);
	commit_init(&com);
	commit_init(&_com);
	for (int i = 0; i < WIDTH; i++) {
		commit_sample_short_crt(r[i]);
	}

	TEST_BEGIN("commitment can be generated and opened") {
		commit_doit(&com, m, &key, r);

		commit_sample_chall_crt(f);
		commit_sample_chall(rho);

		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < NCRT; j++) {
				pcrt_poly_mulmod(s[i][j], r[i][j], f[j], j);
			}
		}

		TEST_ASSERT(commit_open(&com, m, &key, s, f) == 1, end);
	} TEST_END;

	/* With f = 0 and r = 0 the opening equation holds for any message. */
	TEST_ONCE("commitment does not open with a zero challenge") {
		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < NCRT; j++) {
				nmod_poly_zero(_s[i][j]);
			}
		}
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_zero(_f[j]);
		}
		nmod_poly_set(_m, m);
		nmod_poly_set_coeff_ui(_m, 0, (nmod_poly_get_coeff_ui(m, 0) + 1) % MODP);
		TEST_ASSERT(commit_open(&com, _m, &key, _s, _f) == 0, end);
	} TEST_END;

	/* Four coefficients at MODP / 2 and one more chosen so that the sum of
	 * their squares wraps a 64-bit word to 1420468769, far below the bound. */
	TEST_ONCE("norm bound does not wrap on an oversized opening") {
		nmod_poly_zero(_m);
		for (int i = 0; i < 4; i++) {
			nmod_poly_set_coeff_ui(_m, i, MODP / 2);
		}
		nmod_poly_set_coeff_ui(_m, 4, 1785046359);
		TEST_ASSERT(test_norm(_m) == 0, end);
	} TEST_END;

	TEST_BEGIN("commitments are linearly homomorphic") {
		/* Test linearity. */
		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < NCRT; j++) {
				nmod_poly_zero(r[i][j]);
			}
		}
		commit_doit(&_com, rho, &key, r);
		for (int i = 0; i < NCRT; i++) {
			nmod_poly_sub(com.c1[i], com.c1[i], _com.c1[i]);
			nmod_poly_sub(com.c2[i], com.c2[i], _com.c2[i]);
		}
		nmod_poly_sub(m, m, rho);
		TEST_ASSERT(commit_open(&com, m, &key, s, f) == 1, end);
	} TEST_END;

	TEST_BEGIN("sparse product agrees with the full one") {
		nmod_poly_t a, c0, c1;

		nmod_poly_init(a, MODP);
		nmod_poly_init(c0, MODP);
		nmod_poly_init(c1, MODP);
		/* A challenge, and a dense polynomial with a coefficient other than
		 * +-1, which the sparse product has to handle as well. */
		commit_sample_chall(rho);
		nmod_poly_set_coeff_ui(rho, DEGREE - 1, 5);
		nmod_poly_randtest(a, rand, DEGREE);
		commit_poly_mulmod(c0, rho, a);
		commit_poly_mul_sparse(c1, rho, a);
		/* The output may alias an input. */
		commit_poly_mul_sparse(a, rho, a);
		int ok = nmod_poly_equal(c0, c1) && nmod_poly_equal(c0, a);
		nmod_poly_clear(a);
		nmod_poly_clear(c0);
		nmod_poly_clear(c1);
		TEST_ASSERT(ok, end);
	} TEST_END;

  end:
	commit_keyfree(&key);
	commit_free(&com);
	commit_free(&_com);
	nmod_poly_clear(m);
	nmod_poly_clear(_m);
	nmod_poly_clear(rho);
	for (int j = 0; j < NCRT; j++) {
		nmod_poly_clear(f[j]);
		nmod_poly_clear(_f[j]);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_clear(r[i][j]);
			nmod_poly_clear(s[i][j]);
			nmod_poly_clear(_s[i][j]);
		}
	}
}

static void bench(flint_rand_t rand) {
	commitkey_t key;
	commit_t com;
	nmod_poly_t m;
	pcrt_poly_t f, r[WIDTH], s[WIDTH];

	nmod_poly_init(m, MODP);
	nmod_poly_init(f[0], MODP);
	nmod_poly_init(f[1], MODP);
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_init(r[i][j], MODP);
			nmod_poly_init(s[i][j], MODP);
		}
	}

	commit_keyinit(&key);
	commit_keygen(&key);
	commit_init(&com);
	nmod_poly_randtest(m, rand, DEGREE);

	for (int i = 0; i < WIDTH; i++) {
		commit_sample_short_crt(r[i]);
	}

	BENCH_BEGIN("commit_sample") {
		BENCH_ADD(commit_sample_short_crt(r[0]));
	} BENCH_END;

	BENCH_BEGIN("commit_doit") {
		BENCH_ADD(commit_doit(&com, m, &key, r));
	} BENCH_END;

	BENCH_BEGIN("commit_open") {
		commit_sample_chall_crt(f);
		commit_doit(&com, m, &key, r);
		BENCH_ADD(commit_open(&com, m, &key, r, f));
	} BENCH_END;

	commit_keyfree(&key);
	commit_free(&com);
	nmod_poly_clear(m);
	nmod_poly_clear(f[0]);
	nmod_poly_clear(f[1]);
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_clear(r[i][j]);
			nmod_poly_clear(s[i][j]);
		}
	}
}

static void microbench(void) {
	nmod_poly_t alpha, beta, t[2], u[2];

	nmod_poly_init(alpha, MODP);
	nmod_poly_init(beta, MODP);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(t[i], MODP);
		nmod_poly_init(u[i], MODP);
	}

	commit_sample_rand(beta, NULL, DEGREE);
	commit_sample_rand(alpha, NULL, DEGREE);

	BENCH_BEGIN("Polynomial addition") {
		BENCH_ADD(nmod_poly_add(alpha, alpha, beta));
	} BENCH_END;

	BENCH_BEGIN("Polynomial multiplication") {
		BENCH_ADD(commit_poly_mulmod(alpha, alpha, beta));
	} BENCH_END;

	commit_sample_rand_crt(t, NULL);
	commit_sample_rand_crt(u, NULL);

	BENCH_BEGIN("Polynomial mult in CRT form") {
		BENCH_ADD(pcrt_poly_mulmod(t[0], t[0], u[0], 0));
		BENCH_ADD(pcrt_poly_mulmod(t[1], t[1], u[1], 1));
	} BENCH_END;

	nmod_poly_clear(alpha);
	nmod_poly_clear(beta);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(t[i]);
		nmod_poly_clear(u[i]);
	}
}

/* Select which phases to run: "test", "bench", or neither for both. Keeping
 * the benchmarks out of a test run matters in practice, since they dominate
 * the runtime by two orders of magnitude. */
static int phase_selected(int argc, char *argv[], const char *phase) {
	return argc < 2 || strcmp(argv[1], phase) == 0;
}

int main(int argc, char *arv[]) {
	/* Unbuffered, so that a binary that dies still shows which test it
	 * was in: a crash used to take the buffered output with it. */
	setvbuf(stdout, NULL, _IONBF, 0);
	flint_rand_t rand;
	uint64_t buf[2];

	randombytes((unsigned char *)buf, sizeof(buf));
	flint_rand_init(rand);
	flint_rand_set_seed(rand, buf[0], buf[1]);

	commit_setup();

	if (phase_selected(argc, arv, "test")) {
		printf("\n** Tests for lattice-based commitments:\n\n");
		test(rand);
	}

	if (phase_selected(argc, arv, "bench")) {
		printf("\n** Microbenchmarks for polynomial arithmetic:\n\n");
		microbench();

		printf("\n** Benchmarks for lattice-based commitments:\n\n");
		bench(rand);
	}

	commit_finish();
	flint_rand_clear(rand);

	return test_failures() != 0;
}
#endif
