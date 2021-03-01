/**
 * @file
 *
 * Implementation of the lattice-based commitment scheme.
 *
 * @ingroup commit
 */

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

/* Pair of irreducible polynomials for CRT representation. */
static pcrt_poly_t irred;

/* Inverses of the irreducible polynomials for CRT reconstruction. */
static pcrt_poly_t inv;

/**
 * Test if the l2-norm is within bounds (4 * sigma * sqrt(N)).
 *
 * @param[in] r 			- the polynomial to compute the l2-norm.
 * @return the computed norm.
 */
static int test_norm(nmod_poly_t r) {
	// Compute squared norm to save sqrt() and simplify comparison.
	uint64_t norm = commit_norm2_sqr(r);

	// Compute sigma^2 = (11 * v * beta * sqrt(k * N))^2.
	uint64_t sigma_sqr = 11 * NONZERO * BETA;
	sigma_sqr *= sigma_sqr * DEGREE * WIDTH;

	// Compare to (4 * sigma * sqrt(N))^2 = 16 * sigma^2 * N.
	return norm <= (uint64_t)16 * sigma_sqr * DEGREE;
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

// Initialize commitment scheme.
void commit_setup() {
	nmod_poly_init(cyclo_poly, MODP);
	for (int i = 0; i < 2; i++) {
		nmod_poly_init(irred[i], MODP);
		nmod_poly_init(inv[i], MODP);
	}

	// Initialize polynomial as x^N + 1. */
	nmod_poly_set_coeff_ui(cyclo_poly, DEGREE, 1);
	nmod_poly_set_coeff_ui(cyclo_poly, 0, 1);

	// Initialize two factors of the polynomial for CRT representation.
	nmod_poly_set_coeff_ui(irred[0], DEGCRT, 1);
	nmod_poly_set_coeff_ui(irred[0], 0, 3153606543);
	nmod_poly_set_coeff_ui(irred[1], DEGCRT, 1);
	nmod_poly_set_coeff_ui(irred[1], 0, 752843710);

	nmod_poly_invmod(inv[0], irred[0], irred[1]);
	nmod_poly_invmod(inv[1], irred[1], irred[0]);
}

// Finalize commitment scheme.
void commit_finish() {
	for (int i = 0; i < 2; i++) {
		nmod_poly_clear(irred[i]);
		nmod_poly_clear(inv[i]);
	}
	nmod_poly_clear(cyclo_poly);
}

// Return polynomial defining Rp.
nmod_poly_t *commit_poly() {
	return &cyclo_poly;
}

// Return irreducible polynomials defining CRT representation.
nmod_poly_t *commit_irred(int i) {
	return &irred[i];
}

// Recover polynomial from CRT representation.
void pcrt_poly_rec(nmod_poly_t c, pcrt_poly_t a) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);

	nmod_poly_sub(t, a[0], a[1]);
	nmod_poly_mul(t, t, inv[1]);
	nmod_poly_mul(c, t, irred[1]);
	nmod_poly_add(c, c, a[1]);
	nmod_poly_mul(t, irred[0], irred[1]);
	nmod_poly_rem(c, c, t);

	nmod_poly_clear(t);
}

// Compute squared l2-norm.
uint64_t commit_norm2_sqr(nmod_poly_t r) {
	int64_t norm = 0;
	uint32_t coeff;

	/* Compute norm^2. */
	for (int i = 0; i < DEGREE; i++) {
		coeff = nmod_poly_get_coeff_ui(r, i);
		if (coeff > MODP / 2) coeff -= MODP;
		norm += coeff * coeff;
	}
	return norm;
}

// Compute squared l\infty-norm.
uint64_t commit_norm_inf(nmod_poly_t r) {
	int64_t norm = 0;
	uint32_t coeff;

	/* Compute norm_\infty = max absolute coefficient. */
	for (int i = 0; i < DEGREE; i++) {
		coeff = nmod_poly_get_coeff_ui(r, i);
		if (coeff > MODP / 2) coeff = MODP - coeff;
		if (coeff > norm) norm = coeff;
	}
	return norm;
}

// Generate a key pair.
void commit_keygen(commitkey_t *key, flint_rand_t rand) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < 2; k++) {
				nmod_poly_init(key->B1[i][j][k], MODP);
				nmod_poly_zero(key->B1[i][j][k]);
				if (i == j) {
					nmod_poly_set_coeff_ui(key->B1[i][j][k], 0, 1);
				}
			}
		}
	}

	for (int i = 0; i < HEIGHT; i++) {
		for (int j = HEIGHT; j < WIDTH; j++) {
			for (int k = 0; k < 2; k++) {
				nmod_poly_randtest(key->B1[i][j][k], rand, DEGCRT);
			}
		}
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_init(key->b2[i][j], MODP);
			nmod_poly_zero(key->b2[i][j]);
			if (i == HEIGHT) {
				nmod_poly_set_coeff_ui(key->b2[i][j], 0, 1);
			}
			if (i > HEIGHT) {
				nmod_poly_randtest(key->b2[i][j], rand, DEGCRT);
			}
		}
	}
}

// Free a commitment key pair.
void commit_keyfree(commitkey_t *key) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < 2; k++) {
				nmod_poly_clear(key->B1[i][j][k]);
			}
		}
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_clear(key->b2[i][j]);
		}
	}
}

// Sample a short polynomial.
void commit_sample_short(nmod_poly_t r) {
	uint64_t buf;
	uint32_t coeff;

	nmod_poly_zero(r);
	nmod_poly_fit_length(r, DEGREE);
	for (int i = 0; i < DEGREE; i += 32) {
		getrandom(&buf, sizeof(buf), 0);
		for (int j = 0; j < 64; j += 2) {
			if ((buf >> j) & 1) {
				coeff = MODP - ((buf >> (j + 1)) & 1);
			} else {
				coeff = (buf >> (j + 1)) & 1;
			}
			nmod_poly_set_coeff_ui(r, (i+j/2) % DEGREE, coeff);
		}
	}
}

// Sample a short polynomial in CRT representation.
void commit_sample_short_crt(pcrt_poly_t r) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	commit_sample_short(t);
	for (int j = 0; j < 2; j++) {
		nmod_poly_rem(r[j], t, irred[j]);
	}
	nmod_poly_clear(t);
}

// Sample a random polynomial.
void commit_sample_rand(nmod_poly_t r, flint_rand_t rand) {
	nmod_poly_randtest(r, rand, DEGREE);
}

// Sample a random polynomial in CRT representation.
void commit_sample_rand_crt(pcrt_poly_t r, flint_rand_t rand) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	commit_sample_rand(t, rand);
	for (int i = 0; i < 2; i++) {
		nmod_poly_rem(r[i], t, irred[i]);
	}
	nmod_poly_clear(t);
}

// Sample a challenge.
void commit_sample_chall(nmod_poly_t f) {
	nmod_poly_zero(f);
	nmod_poly_t c[2];
	uint32_t buf;

	for (int i = 0; i < 2; i++) {
		nmod_poly_init(c[i], MODP);
		nmod_poly_fit_length(c[i], DEGREE);
		for (int j = 0; j < NONZERO; j++) {
			getrandom(&buf, sizeof(buf), 0);
			buf = buf % DEGREE;
			while (nmod_poly_get_coeff_ui(c[i], buf) != 0) {
				getrandom(&buf, sizeof(buf), 0);
				buf = buf % DEGREE;
			}
			nmod_poly_set_coeff_ui(c[i], buf, 1);
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
	nmod_poly_rem(f[0], t, irred[0]);
	nmod_poly_rem(f[1], t, irred[1]);
}

// Sample a polynomial according to a Gaussian distribution.
void commit_sample_gauss(nmod_poly_t r) {
	int32_t coeff;
	for (int i = 0; i < DEGREE; i ++) {
		coeff = discrete_gaussian(0.0);
		if (coeff < 0) coeff += MODP;
		nmod_poly_set_coeff_ui(r, i, coeff);
	}
}

// Sample a polynomial according to a Gaussian distribution in CRT rep.
void commit_sample_gauss_crt(nmod_poly_t r[2]) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	commit_sample_gauss(t);
	nmod_poly_rem(r[0], t, irred[0]);
	nmod_poly_rem(r[1], t, irred[1]);

	nmod_poly_clear(t);
}

// Commit to a message.
void commit_doit(commit_t *com, nmod_poly_t m, commitkey_t *key, pcrt_poly_t r[WIDTH]) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	for (int i = 0; i < 2; i++) {
		nmod_poly_init(com->c1[i], MODP);
		nmod_poly_init(com->c2[i], MODP);
		nmod_poly_zero(com->c1[i]);
		nmod_poly_zero(com->c2[i]);
	}

	// Compute B = [ B1 b2 ]^t * r_m.
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < 2; k++) {
				nmod_poly_mulmod(t, key->B1[i][j][k], r[j][k], irred[k]);
				nmod_poly_add(com->c1[k], com->c1[k], t);
				nmod_poly_mulmod(t, key->b2[j][k], r[j][k], irred[k]);
				nmod_poly_add(com->c2[k], com->c2[k], t);
			}
		}
	}

	// Convert m to CRT representation and accumulate.
	for (int i = 0; i < 2; i++) {
		nmod_poly_rem(t, m, irred[i]);
		nmod_poly_add(com->c2[i], com->c2[i], t);
	}

	nmod_poly_clear(t);
}

// Open a commitment on a message, randomness, factor.
int commit_open(commit_t *com, nmod_poly_t m, commitkey_t *key,  pcrt_poly_t r[WIDTH], pcrt_poly_t f) {
	nmod_poly_t t;
	pcrt_poly_t c1, c2, _c1, _c2;
	int result = 0;

	nmod_poly_init(t, MODP);
	for (int i = 0; i < 2; i++) {
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
			for (int k = 0; k < 2; k++) {
				nmod_poly_mulmod(t, key->B1[i][j][k], r[j][k], irred[k]);
				nmod_poly_add(c1[k], c1[k], t);
				nmod_poly_mulmod(t, key->b2[j][k], r[j][k], irred[k]);
				nmod_poly_add(c2[k], c2[k], t);
			}
		}
	}

	// Convert m to CRT representation before multiplication.
	for (int i = 0; i < 2; i++) {
		nmod_poly_rem(t, m, irred[i]);
		nmod_poly_mulmod(t, t, f[i], irred[i]);
		nmod_poly_add(c2[i], c2[i], t);
	}

	for (int i = 0; i < 2; i++) {
		nmod_poly_mulmod(_c1[i], com->c1[i], f[i], irred[i]);
		nmod_poly_mulmod(_c2[i], com->c2[i], f[i], irred[i]);
	}

	pcrt_poly_rec(t, r[0]);
	if (test_norm(t)) {
		pcrt_poly_rec(t, r[1]);
		if (test_norm(t)) {
			pcrt_poly_rec(t, r[2]);
			if (test_norm(t)) {
				if (nmod_poly_equal(_c1[0], c1[0]) && nmod_poly_equal(_c2[0], c2[0])) {
					if (nmod_poly_equal(_c1[1], c1[1]) && nmod_poly_equal(_c2[1], c2[1])) {
						result = 1;
					}
				}
			}
		}
	}

	nmod_poly_clear(t);
	for (int i = 0; i < 2; i++) {
		nmod_poly_clear(c1[i]);
		nmod_poly_clear(c2[i]);
	}
	return result;
}

// Free a commitment.
void commit_free(commit_t *com) {
	for (int i = 0; i < 2; i++) {
		nmod_poly_clear(com->c1[i]);
		nmod_poly_clear(com->c2[i]);
	}
}

#ifdef MAIN
// Tests and benchmarks below.
static void test(flint_rand_t rand) {
	commitkey_t key;
	commit_t com, _com;
	nmod_poly_t m, rho;
	pcrt_poly_t r[WIDTH], s[WIDTH], f;

	nmod_poly_init(m, MODP);
	nmod_poly_init(rho, MODP);
	nmod_poly_init(f[0], MODP);
	nmod_poly_init(f[1], MODP);
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_init(r[i][j], MODP);
			nmod_poly_init(s[i][j], MODP);
		}
	}

	/* Generate a random message. */
	nmod_poly_randtest(m, rand, DEGREE);

	/* Generate commitment key. */
	commit_keygen(&key, rand);
	for (int i = 0; i < WIDTH; i++) {
		commit_sample_short_crt(r[i]);
	}

	TEST_BEGIN("commitment can be generated and opened") {

		commit_doit(&com, m, &key, r);

		commit_sample_chall_crt(f);
		commit_sample_chall(rho);

		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < 2; j++) {
				nmod_poly_mulmod(s[i][j], r[i][j], f[j], irred[j]);
			}
		}

		TEST_ASSERT(commit_open(&com, m, &key, s, f) == 1, end);
	} TEST_END;

	TEST_BEGIN("commitments are linearly homomorphic") {
		/* Test linearity. */
		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < 2; j++) {
				nmod_poly_zero(r[i][j]);
			}
		}
		commit_doit(&_com, rho, &key, r);
		for (int i = 0; i < 2; i++) {
			nmod_poly_sub(com.c1[i], com.c1[i], _com.c1[i]);
			nmod_poly_sub(com.c2[i], com.c2[i], _com.c2[i]);
		}
		nmod_poly_sub(m, m, rho);
		TEST_ASSERT(commit_open(&com, m, &key, s, f) == 1, end);
	} TEST_END;

end:
	commit_keyfree(&key);
	commit_free(&com);
	nmod_poly_clear(m);
	nmod_poly_clear(rho);
	nmod_poly_clear(f[0]);
	nmod_poly_clear(f[1]);
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_clear(r[i][j]);
			nmod_poly_clear(s[i][j]);
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
		for (int j = 0; j < 2; j++) {
			nmod_poly_init(r[i][j], MODP);
			nmod_poly_init(s[i][j], MODP);
		}
	}

	commit_keygen(&key, rand);
	nmod_poly_randtest(m, rand, DEGREE);

	for (int i = 0; i < WIDTH; i++) {
		commit_sample_short_crt(r[i]);
	}

	BENCH_BEGIN("commit_sample") {
		BENCH_ADD(commit_sample_short_crt(r[0]));
	} BENCH_END;

	BENCH_BEGIN("commit_doit") {
		BENCH_ADD(commit_doit(&com, m, &key, r));
		commit_free(&com);
	} BENCH_END;

	BENCH_BEGIN("commit_open") {
		commit_sample_chall_crt(f);
		commit_doit(&com, m, &key, r);
		BENCH_ADD(commit_open(&com, m, &key, s, f));
	} BENCH_END;

	commit_keyfree(&key);
	commit_free(&com);
	nmod_poly_clear(m);
	nmod_poly_clear(f[0]);
	nmod_poly_clear(f[1]);
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_clear(r[i][j]);
			nmod_poly_clear(s[i][j]);
		}
	}
}

int main(int argc, char *arv[]) {
	flint_rand_t rand;

	flint_randinit(rand);
	commit_setup();

	printf("\n** Tests for lattice-based commitments:\n\n");
	test(rand);

	printf("\n** Benchmarks for lattice-based commitments:\n\n");
	bench(rand);

	commit_finish();
}
#endif
