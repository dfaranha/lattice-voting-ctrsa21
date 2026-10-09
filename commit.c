/**
 * @file
 *
 * Implementation of the lattice-based commitment scheme.
 *
 * @ingroup commit
 */

#include <flint/nmod.h>

#include "param.h"
#include "commit.h"
#include "test.h"
#include "bench.h"
#include "gaussian.h"
#include "fastrandombytes.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/* The first square root of -1. */
#define	P0		3153606543
/* The second square root of -1. */
#define P1		752843710

/* Polynomial defining the cyclotomic ring. */
static nmod_poly_t cyclo_poly;

/* Pair of irreducible polynomials that (x^N + 1) splits into. They are kept so
 * that commit_irred still describes the ring, but this build does not use them
 * for arithmetic. */
static nmod_poly_t irred[NCRT];

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
	 * accepted only openings 2.5 times shorter than extraction yields. */
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

// Initialize commitment scheme.
void commit_setup() {
	nmod_poly_init(cyclo_poly, MODP);
	nmod_poly_init(mul_tmp, MODP);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(irred[i], MODP);
	}

	// Initialize polynomial as x^N + 1. */
	nmod_poly_set_coeff_ui(cyclo_poly, DEGREE, 1);
	nmod_poly_set_coeff_ui(cyclo_poly, 0, 1);

	// Initialize two factors of the polynomial for CRT representation.
	nmod_poly_set_coeff_ui(irred[0], DEGCRT, 1);
	nmod_poly_set_coeff_ui(irred[0], 0, 3153606543);
	nmod_poly_set_coeff_ui(irred[1], DEGCRT, 1);
	nmod_poly_set_coeff_ui(irred[1], 0, 752843710);

}

// Finalize commitment scheme.
void commit_finish() {
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(irred[i]);
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


// Multiply two polynomials in Rp.
void commit_poly_mulmod(nmod_poly_t c, const nmod_poly_t a,
		const nmod_poly_t b) {
	nmod_poly_mul(mul_tmp, a, b);
	fold_mod(c, mul_tmp, DEGREE, 1);
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
		 * coefficient below it can still square past 2^64 and wrap, which
		 * would let an oversized response pass this check. Dividing decides
		 * the same question without overflowing. Wrapping needs a centred
		 * coefficient of 2^32, so at the current MODP this is unreachable,
		 * with about one bit to spare. */
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
			nmod_poly_init(key->B1[i][j], MODP);
		}
	}
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(key->b2[i], MODP);
	}
}

// Generate a key pair.
void commit_keygen(commitkey_t *key) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			nmod_poly_zero(key->B1[i][j]);
			if (i == j) {
				nmod_poly_set_coeff_ui(key->B1[i][j], 0, 1);
			}
		}
	}

	for (int i = 0; i < HEIGHT; i++) {
		for (int j = HEIGHT; j < WIDTH; j++) {
			/* A uniform ring element spans the whole degree here; the CRT
			 * build drew DEGCRT coefficients for each of its components. */
			commit_sample_rand(key->B1[i][j], NULL, DEGREE);
		}
	}
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_zero(key->b2[i]);
		if (i == HEIGHT) {
			nmod_poly_set_coeff_ui(key->b2[i], 0, 1);
		}
		if (i > HEIGHT) {
			commit_sample_rand(key->b2[i], NULL, DEGREE);
		}
	}
}

// Free a commitment key pair.
void commit_keyfree(commitkey_t *key) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			nmod_poly_clear(key->B1[i][j]);
		}
	}
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_clear(key->b2[i]);
	}
}

// Sample a short polynomial.
void commit_sample_short(nmod_poly_t r) {
	uint64_t buf;
	/* This holds MODP - 1 + d for d in {0, 1, 2}, so it has to be as wide as
	 * the modulus. As a uint32_t it is correct at the current MODP but only
	 * just: the moment MODP passes 2^32 it truncates, the randomness stops
	 * being ternary, and rejection sampling then never accepts, so the prover
	 * hangs rather than failing. */
	ulong coeff;
	int i, j, s;

	nmod_poly_zero(r);
	nmod_poly_fit_length(r, DEGREE);
	i = 0;
	s = 8 * sizeof(buf);
	j = s;

	do {
		if (j == s) {
			getrandom(&buf, sizeof(buf), 0);
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
				getrandom(&buf, sizeof(buf), 0);
				pos = buf % DEGREE;
			} while (nmod_poly_get_coeff_ui(c[i], pos) != 0);
			nmod_poly_set_coeff_ui(c[i], pos, (buf >> 31) ? MODP - 1 : 1);
		}
	}
	nmod_poly_sub(f, c[0], c[1]);

	nmod_poly_clear(c[0]);
	nmod_poly_clear(c[1]);
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


// Initialise a commitment.
void commit_init(commit_t *com) {
	nmod_poly_init(com->c1, MODP);
	nmod_poly_init(com->c2, MODP);
}

// Commit to a message.
void commit_doit(commit_t *com, nmod_poly_t m, commitkey_t *key,
		nmod_poly_t r[WIDTH]) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	nmod_poly_zero(com->c1);
	nmod_poly_zero(com->c2);

	// Compute B = [ B1 b2 ]^t * r_m.
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			commit_poly_mulmod(t, key->B1[i][j], r[j]);
			nmod_poly_add(com->c1, com->c1, t);
			if (i == 0) {
				commit_poly_mulmod(t, key->b2[j], r[j]);
				nmod_poly_add(com->c2, com->c2, t);
			}
		}
	}

	nmod_poly_add(com->c2, com->c2, m);

	nmod_poly_clear(t);
}

// Open a commitment on a message, randomness, factor.
int commit_open(commit_t *com, nmod_poly_t m, commitkey_t *key,
		nmod_poly_t r[WIDTH], nmod_poly_t f) {
	nmod_poly_t t, c1, c2, _c1, _c2;
	int result = 0;

	nmod_poly_init(t, MODP);
	nmod_poly_init(c1, MODP);
	nmod_poly_init(c2, MODP);
	nmod_poly_init(_c1, MODP);
	nmod_poly_init(_c2, MODP);
	nmod_poly_zero(c1);
	nmod_poly_zero(c2);

	// Compute B = [ B1 b2 ]^t * r_m.
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			commit_poly_mulmod(t, key->B1[i][j], r[j]);
			nmod_poly_add(c1, c1, t);
			if (i == 0) {
				commit_poly_mulmod(t, key->b2[j], r[j]);
				nmod_poly_add(c2, c2, t);
			}
		}
	}

	commit_poly_mulmod(t, m, f);
	nmod_poly_add(c2, c2, t);

	commit_poly_mulmod(_c1, com->c1, f);
	commit_poly_mulmod(_c2, com->c2, f);

	/* The challenge and randomness are already ring elements, so there is
	 * nothing to reconstruct before checking them. Every component of r is
	 * bounded, not only the first three, which is all this used to check
	 * whatever the WIDTH. */
	result = test_chall_diff(f);
	for (int i = 0; i < WIDTH; i++) {
		result &= test_norm(r[i]);
	}
	result &= nmod_poly_equal(_c1, c1);
	result &= nmod_poly_equal(_c2, c2);

	nmod_poly_clear(t);
	nmod_poly_clear(c1);
	nmod_poly_clear(c2);
	nmod_poly_clear(_c1);
	nmod_poly_clear(_c2);
	return result;
}

// Free a commitment.
void commit_free(commit_t *com) {
	nmod_poly_clear(com->c1);
	nmod_poly_clear(com->c2);
}

#ifdef MAIN
// Tests and benchmarks below.
static void test(flint_rand_t rand) {
	commitkey_t key;
	commit_t com, _com;
	nmod_poly_t m, _m, rho;
	nmod_poly_t r[WIDTH], s[WIDTH], _s[WIDTH], f, _f;

	nmod_poly_init(m, MODP);
	nmod_poly_init(_m, MODP);
	nmod_poly_init(rho, MODP);
	nmod_poly_init(f, MODP);
	nmod_poly_init(_f, MODP);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(r[i], MODP);
		nmod_poly_init(s[i], MODP);
		nmod_poly_init(_s[i], MODP);
	}

	/* Generate a random message. */
	nmod_poly_randtest(m, rand, DEGREE);

	/* Generate commitment key. */
	commit_keyinit(&key);
	commit_keygen(&key);
	commit_init(&com);
	commit_init(&_com);
	for (int i = 0; i < WIDTH; i++) {
		commit_sample_short(r[i]);
	}

	TEST_BEGIN("commitment can be generated and opened") {
		commit_doit(&com, m, &key, r);

		commit_sample_chall(f);
		commit_sample_chall(rho);

		for (int i = 0; i < WIDTH; i++) {
			commit_poly_mulmod(s[i], r[i], f);
		}

		TEST_ASSERT(commit_open(&com, m, &key, s, f) == 1, end);
	} TEST_END;

	/* With f = 0 and r = 0 the opening equation holds for any message. */
	TEST_ONCE("commitment does not open with a zero challenge") {
		for (int i = 0; i < WIDTH; i++) {
			nmod_poly_zero(_s[i]);
		}
		nmod_poly_zero(_f);
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
			nmod_poly_zero(r[i]);
		}
		commit_doit(&_com, rho, &key, r);
		nmod_poly_sub(com.c1, com.c1, _com.c1);
		nmod_poly_sub(com.c2, com.c2, _com.c2);
		nmod_poly_sub(m, m, rho);
		TEST_ASSERT(commit_open(&com, m, &key, s, f) == 1, end);
	} TEST_END;

  end:
	commit_keyfree(&key);
	commit_free(&com);
	commit_free(&_com);
	nmod_poly_clear(m);
	nmod_poly_clear(_m);
	nmod_poly_clear(rho);
	nmod_poly_clear(f);
	nmod_poly_clear(_f);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_clear(r[i]);
		nmod_poly_clear(s[i]);
		nmod_poly_clear(_s[i]);
	}
}

static void bench(flint_rand_t rand) {
	commitkey_t key;
	commit_t com;
	nmod_poly_t m;
	nmod_poly_t f, r[WIDTH], s[WIDTH];

	nmod_poly_init(m, MODP);
	nmod_poly_init(f, MODP);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(r[i], MODP);
		nmod_poly_init(s[i], MODP);
	}

	commit_keyinit(&key);
	commit_keygen(&key);
	commit_init(&com);
	nmod_poly_randtest(m, rand, DEGREE);

	for (int i = 0; i < WIDTH; i++) {
		commit_sample_short(r[i]);
	}

	BENCH_BEGIN("commit_sample") {
		BENCH_ADD(commit_sample_short(r[0]));
	} BENCH_END;

	BENCH_BEGIN("commit_doit") {
		BENCH_ADD(commit_doit(&com, m, &key, r));
	} BENCH_END;

	BENCH_BEGIN("commit_open") {
		commit_sample_chall(f);
		commit_doit(&com, m, &key, r);
		BENCH_ADD(commit_open(&com, m, &key, r, f));
	} BENCH_END;

	commit_keyfree(&key);
	commit_free(&com);
	nmod_poly_clear(m);
	nmod_poly_clear(f);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_clear(r[i]);
		nmod_poly_clear(s[i]);
	}
}

static void microbench(void) {
	nmod_poly_t alpha, beta;

	nmod_poly_init(alpha, MODP);
	nmod_poly_init(beta, MODP);

	commit_sample_rand(beta, NULL, DEGREE);
	commit_sample_rand(alpha, NULL, DEGREE);

	BENCH_BEGIN("Polynomial addition") {
		BENCH_ADD(nmod_poly_add(alpha, alpha, beta));
	} BENCH_END;

	/* This build multiplies directly in Rp; there is no CRT form to compare
	 * against here. See the crt branch for the split representation. */
	BENCH_BEGIN("Polynomial multiplication") {
		BENCH_ADD(commit_poly_mulmod(alpha, alpha, beta));
	} BENCH_END;

	nmod_poly_clear(alpha);
	nmod_poly_clear(beta);
}

/* Select which phases to run: "test", "bench", or neither for both. Keeping
 * the benchmarks out of a test run matters in practice, since they dominate
 * the runtime by two orders of magnitude. */
static int phase_selected(int argc, char *argv[], const char *phase) {
	return argc < 2 || strcmp(argv[1], phase) == 0;
}

int main(int argc, char *arv[]) {
	flint_rand_t rand;
	uint64_t buf[2];

	getrandom(buf, sizeof(buf), GRND_RANDOM);
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
