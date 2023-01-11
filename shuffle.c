#include <math.h>
#include <stdlib.h>

#include "param.h"
#include "commit.h"
#include "test.h"
#include "bench.h"
#include "assert.h"
#include "fastrandombytes.h"
#include "sha.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#define MSGS        25

void shuffle_hash(nmod_poly_t beta, commit_t c[MSGS], commit_t d[MSGS]) {
	flint_rand_t rand;
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	uint64_t seed0, seed1;

	SHA256Reset(&sha);

	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < 2; j++) {
			SHA256Input(&sha, (const uint8_t*)c[i].c1[j]->coeffs, c[i].c1[j]->alloc * sizeof(uint64_t));
			SHA256Input(&sha, (const uint8_t*)c[i].c2[j]->coeffs, c[i].c2[j]->alloc * sizeof(uint64_t));
			SHA256Input(&sha, (const uint8_t*)d[i].c1[j]->coeffs, d[i].c1[j]->alloc * sizeof(uint64_t));
			SHA256Input(&sha, (const uint8_t*)d[i].c2[j]->coeffs, d[i].c2[j]->alloc * sizeof(uint64_t));
		}
	}

	SHA256Result(&sha, hash);

	flint_randinit(rand);
	memcpy(&seed0, hash, sizeof(uint64_t));
	memcpy(&seed1, hash + SHA256HashSize/2, sizeof(uint64_t));
	flint_randseed(rand, seed0, seed1);
	commit_sample_rand(beta, rand);
	flint_randclear(rand);
}

void lin_hash(nmod_poly_t d[2], commitkey_t *key, commit_t x, commit_t y,
		nmod_poly_t alpha, nmod_poly_t beta, nmod_poly_t u[2],
		nmod_poly_t t[2], nmod_poly_t _t[2]) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	uint32_t buf;

	SHA256Reset(&sha);

	/* Hash public key. */
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < 2; k++) {
				SHA256Input(&sha, (const uint8_t*)key->B1[i][j][k]->coeffs, key->B1[i][j][k]->alloc * sizeof(uint64_t));
				if (i == 0) {
					SHA256Input(&sha, (const uint8_t*)key->b2[j][k]->coeffs, key->b2[j][k]->alloc * sizeof(uint64_t));
				}
			}
		}
	}

	/* Hash alpha, beta from linear relation. */
	SHA256Input(&sha, (const uint8_t*)alpha->coeffs, alpha->alloc * sizeof(uint64_t));
	SHA256Input(&sha, (const uint8_t*)beta->coeffs, beta->alloc * sizeof(uint64_t));

	/* Hash [x], [x'], t, t' in CRT representation. */
	for (int i = 0; i < 2; i++) {
		SHA256Input(&sha, (const uint8_t*)x.c1[i]->coeffs, x.c1[i]->alloc * sizeof(uint64_t));
		SHA256Input(&sha, (const uint8_t*)x.c2[i]->coeffs, x.c2[i]->alloc * sizeof(uint64_t));
		SHA256Input(&sha, (const uint8_t*)y.c1[i]->coeffs, y.c1[i]->alloc * sizeof(uint64_t));
		SHA256Input(&sha, (const uint8_t*)y.c2[i]->coeffs, y.c2[i]->alloc * sizeof(uint64_t));
		SHA256Input(&sha, (const uint8_t*)u[i]->coeffs, u[i]->alloc * sizeof(uint64_t));
		SHA256Input(&sha, (const uint8_t*)t[i]->coeffs, t[i]->alloc * sizeof(uint64_t));
		SHA256Input(&sha, (const uint8_t*)_t[i]->coeffs, _t[i]->alloc * sizeof(uint64_t));
	}

	SHA256Result(&sha, hash);

	/* Sample challenge from RNG seeded with hash. */
	fastrandombytes_setseed(hash);
	for (int i = 0; i < 2; i++) {
		nmod_poly_fit_length(d[i], DEGREE);
		for (int j = 0; j < NONZERO; j++) {
			fastrandombytes((unsigned char *)&buf, sizeof(buf));
			buf = buf % DEGREE;
			while (nmod_poly_get_coeff_ui(d[i], buf) != 0) {
				fastrandombytes((unsigned char *)&buf, sizeof(buf));
				buf = buf % DEGREE;
			}
			nmod_poly_set_coeff_ui(d[i], buf, 1);
		}
	}
	nmod_poly_sub(d[1], d[0], d[1]);
	nmod_poly_rem(d[0], d[1], *commit_irred(0));
	nmod_poly_rem(d[1], d[1], *commit_irred(1));
	/* TODO: rejection sampling. */
}

static void lin_prover(nmod_poly_t y[WIDTH][2], nmod_poly_t _y[WIDTH][2],
		nmod_poly_t t[2], nmod_poly_t _t[2], nmod_poly_t u[2],
		commit_t x, commit_t _x, commitkey_t *key, nmod_poly_t alpha,
		nmod_poly_t beta, nmod_poly_t r[WIDTH][2], nmod_poly_t _r[WIDTH][2],
		int l) {
	nmod_poly_t tmp, d[2];

	nmod_poly_init(tmp, MODP);
	for (int i = 0; i < 2; i++) {
		nmod_poly_init(d[i], MODP);
		nmod_poly_zero(t[i]);
		nmod_poly_zero(_t[i]);
		nmod_poly_zero(u[i]);
	}

	for (int i = 0; i < WIDTH; i++) {
		commit_sample_gauss_crt(y[i]);
		commit_sample_gauss_crt(_y[i]);
	}
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < 2; k++) {
				nmod_poly_mulmod(tmp, key->B1[i][j][k], y[j][k], *commit_irred(k));
				nmod_poly_add(t[k], t[k], tmp);
				nmod_poly_mulmod(tmp, key->B1[i][j][k], _y[j][k], *commit_irred(k));
				nmod_poly_add(_t[k], _t[k], tmp);
			}
		}
	}

	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_mulmod(tmp, key->b2[i][j], y[i][j], *commit_irred(j));
			nmod_poly_mulmod(tmp, tmp, alpha, *commit_irred(j));
			nmod_poly_add(u[j], u[j], tmp);
			nmod_poly_mulmod(tmp, key->b2[i][j], _y[i][j], *commit_irred(j));
			nmod_poly_sub(u[j], u[j], tmp);
		}
	}

	/* Sample challenge. */
	lin_hash(d, key, x, _x, alpha, beta, u, t, _t);

	/* Prover */
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_mulmod(tmp, d[j], r[i][j], *commit_irred(j));
			nmod_poly_add(y[i][j], y[i][j], tmp);
			nmod_poly_mulmod(tmp, d[j], _r[i][j], *commit_irred(j));
			nmod_poly_add(_y[i][j], _y[i][j], tmp);
		}
	}

	nmod_poly_clear(tmp);
	for (int i = 0; i < 2; i++) {
		nmod_poly_clear(d[i]);
	}
}

static int lin_verifier(nmod_poly_t y[WIDTH][2], nmod_poly_t _y[WIDTH][2],
		nmod_poly_t t[2], nmod_poly_t _t[2], nmod_poly_t u[2],
		commit_t com, commit_t x, commitkey_t *key,
		nmod_poly_t alpha, nmod_poly_t beta, int l) {
	nmod_poly_t tmp, _d[2], v[2], _v[2], z[WIDTH], _z[WIDTH];
	int result = 1;

	nmod_poly_init(tmp, MODP);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(z[i], MODP);
		nmod_poly_init(_z[i], MODP);
	}
	for (int i = 0; i < 2; i++) {
		nmod_poly_init(_d[i], MODP);
		nmod_poly_init(v[i], MODP);
		nmod_poly_init(_v[i], MODP);
		nmod_poly_zero(v[i]);
		nmod_poly_zero(_v[i]);
	}

	/* Sample challenge. */
	lin_hash(_d, key, com, x, alpha, beta, u, t, _t);

	/* Verifier checks norm, reconstruct from CRT representation. */
	for (int i = 0; i < WIDTH; i++) {
		pcrt_poly_rec(z[i], y[i]);
		pcrt_poly_rec(_z[i], _y[i]);
		assert(commit_norm2_sqr(z[i]) <= (uint64_t)4 * DEGREE * SIGMA_C * SIGMA_C);
		assert(commit_norm2_sqr(_z[i]) <= (uint64_t)4 * DEGREE * SIGMA_C * SIGMA_C);
	}
	/* Verifier computes B1z and B1z'. */
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < 2; k++) {
				nmod_poly_mulmod(tmp, key->B1[i][j][k], y[j][k], *commit_irred(k));
				nmod_poly_add(v[k], v[k], tmp);
				nmod_poly_mulmod(tmp, key->B1[i][j][k], _y[j][k], *commit_irred(k));
				nmod_poly_add(_v[k], _v[k], tmp);
			}
		}
	}
	/* Verifier checks that B_1z = t + dc1, B_1z' = t' + dc1'. */
	for (int j = 0; j < 2; j++) {
		nmod_poly_mulmod(tmp, _d[j], com.c1[j], *commit_irred(j));
		nmod_poly_add(t[j], t[j], tmp);
		nmod_poly_mulmod(tmp, _d[j], x.c1[j], *commit_irred(j));
		nmod_poly_add(_t[j], _t[j], tmp);
		result &= nmod_poly_equal(t[j], v[j]);
		result &= nmod_poly_equal(_t[j], _v[j]);
	}

	if (l == 0) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_mulmod(t[j], alpha, com.c2[j], *commit_irred(j));
			nmod_poly_add(t[j], t[j], beta);
			nmod_poly_sub(t[j], t[j], x.c2[j]);
			nmod_poly_mulmod(t[j], t[j], _d[j], *commit_irred(j));
			nmod_poly_add(t[j], t[j], u[j]);
			nmod_poly_rem(t[j], t[j], *commit_irred(j));
		}
	}
	if (l > 0 && l < MSGS - 1) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_mulmod(t[j], alpha, com.c2[j], *commit_irred(j));
			nmod_poly_add(t[j], t[j], beta);
			nmod_poly_sub(t[j], t[j], x.c2[j]);
			nmod_poly_mulmod(t[j], t[j], _d[j], *commit_irred(j));
			nmod_poly_add(t[j], t[j], u[j]);
		}
	}
	if (l == MSGS - 1) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_mulmod(t[j], alpha, com.c2[j], *commit_irred(j));
			if (MSGS & 1) {
				nmod_poly_sub(t[j], t[j], beta);
			} else {
				nmod_poly_add(t[j], t[j], beta);
			}
			nmod_poly_sub(t[j], t[j], x.c2[j]);
			nmod_poly_mulmod(t[j], t[j], _d[j], *commit_irred(j));
			nmod_poly_add(t[j], t[j], u[j]);
		}
	}

	nmod_poly_zero(v[0]);
	nmod_poly_zero(v[1]);
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_mulmod(tmp, key->b2[i][j], y[i][j], *commit_irred(j));
			nmod_poly_mulmod(tmp, alpha, tmp, *commit_irred(j));
			nmod_poly_add(v[j], v[j], tmp);
			nmod_poly_mulmod(tmp, key->b2[i][j], _y[i][j], *commit_irred(j));
			nmod_poly_sub(v[j], v[j], tmp);
		}
	}
	for (int j = 0; j < 2; j++) {
		result &= nmod_poly_equal(t[j], v[j]);
	}

	nmod_poly_clear(tmp);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_clear(z[i]);
		nmod_poly_clear(_z[i]);
	}
	for (int i = 0; i < 2; i++) {
		nmod_poly_clear(_d[i]);
		nmod_poly_clear(v[i]);
		nmod_poly_clear(_v[i]);
	}
	return result;
}

static void shuffle_prover(nmod_poly_t y[MSGS][WIDTH][2], nmod_poly_t _y[MSGS][WIDTH][2],
		nmod_poly_t t[MSGS][2], nmod_poly_t _t[MSGS][2], nmod_poly_t u[MSGS][2],
		commit_t d[MSGS], nmod_poly_t s[MSGS], commit_t com[MSGS],
		nmod_poly_t m[MSGS], nmod_poly_t _m[MSGS], nmod_poly_t r[MSGS][WIDTH][2],
		nmod_poly_t rho, commitkey_t *key, flint_rand_t rng) {
	nmod_poly_t beta, t0, t1, theta[MSGS], _r[MSGS][WIDTH][2];

	nmod_poly_init(t0, MODP);
	nmod_poly_init(t1, MODP);
	nmod_poly_init(beta, MODP);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(theta[i], MODP);
		for (int k = 0; k < 2; k++) {
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_init(_r[i][j][k], MODP);
			}
		}
	}

	/* Prover shifts the messages by rho. */
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_sub(m[i], m[i], rho);
		nmod_poly_sub(_m[i], _m[i], rho);
	}

	/* Prover samples theta_i and computes commitments D_i. */
	commit_sample_rand(theta[0], rng);
	nmod_poly_mulmod(t0, theta[0], _m[0], *commit_poly());
	for (int j = 0; j < WIDTH; j++) {
		commit_sample_short_crt(_r[0][j]);
	}
	commit_doit(&d[0], t0, key, _r[0]);
	for (int i = 1; i < MSGS - 1; i++) {
		commit_sample_rand(theta[i], rng);
		nmod_poly_mulmod(t0, theta[i - 1], m[i], *commit_poly());
		nmod_poly_mulmod(t1, theta[i], _m[i], *commit_poly());
		nmod_poly_add(t0, t0, t1);
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short_crt(_r[i][j]);
		}
		commit_doit(&d[i], t0, key, _r[i]);
	}
	nmod_poly_mulmod(t0, theta[MSGS - 2], m[MSGS - 1], *commit_poly());
	for (int j = 0; j < WIDTH; j++) {
		commit_sample_short_crt(_r[MSGS - 1][j]);
	}
	commit_doit(&d[MSGS - 1], t0, key, _r[MSGS - 1]);

	shuffle_hash(beta, com, d);
	nmod_poly_mulmod(s[0], theta[0], _m[0], *commit_poly());
	nmod_poly_mulmod(t0, beta, m[0], *commit_poly());
	nmod_poly_sub(s[0], s[0], t0);
	nmod_poly_invmod(t0, _m[0], *commit_poly());
	nmod_poly_mulmod(s[0], s[0], t0, *commit_poly());
	for (int i = 1; i < MSGS - 1; i++) {
		nmod_poly_mulmod(s[i], theta[i - 1], m[i], *commit_poly());
		nmod_poly_mulmod(t0, theta[i], _m[i], *commit_poly());
		nmod_poly_add(s[i], s[i], t0);
		nmod_poly_mulmod(t0, s[i - 1], m[i], *commit_poly());
		nmod_poly_sub(s[i], s[i], t0);
		nmod_poly_invmod(t0, _m[i], *commit_poly());
		nmod_poly_mulmod(s[i], s[i], t0, *commit_poly());
	}

	for (int l = 0; l < MSGS; l++) {
		if (l < MSGS - 1) {
			nmod_poly_mulmod(t0, s[l], _m[l], *commit_poly());
		} else {
			nmod_poly_mulmod(t0, beta, _m[l], *commit_poly());
		}

		if (l == 0) {
			lin_prover(y[l], _y[l], t[l], _t[l], u[l], com[l], d[l], key, beta, t0, r[l], _r[l], l);
		} else {
			lin_prover(y[l], _y[l], t[l], _t[l], u[l], com[l], d[l], key, s[l - 1], t0, r[l], _r[l], l);
		}
	}

	nmod_poly_clear(t0);
	nmod_poly_clear(t1);
	nmod_poly_clear(beta);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_clear(theta[i]);
		for (int k = 0; k < 2; k++) {
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_clear(_r[i][j][k]);
			}
		}
	}
}

static int shuffle_verifier(nmod_poly_t y[MSGS][WIDTH][2], nmod_poly_t _y[MSGS][WIDTH][2],
		nmod_poly_t t[MSGS][2], nmod_poly_t _t[MSGS][2], nmod_poly_t u[MSGS][2],
		commit_t d[MSGS], nmod_poly_t s[MSGS], commit_t com[MSGS],
		nmod_poly_t _m[MSGS], commitkey_t *key) {
	int result = 1;
	nmod_poly_t beta, t0;

	nmod_poly_init(t0, MODP);
	nmod_poly_init(beta, MODP);

	shuffle_hash(beta, com, d);
	/* Now verify each \Prod_LIN instance, one for each commitment. */
	for (int l = 0; l < MSGS; l++) {
		if (l < MSGS - 1) {
			nmod_poly_mulmod(t0, s[l], _m[l], *commit_poly());
		} else {
			nmod_poly_mulmod(t0, beta, _m[l], *commit_poly());
		}

		if (l == 0) {
			result &= lin_verifier(y[l], _y[l], t[l], _t[l], u[l], com[l], d[l], key, beta, t0, l);
		} else {
			result &= lin_verifier(y[l], _y[l], t[l], _t[l], u[l], com[l], d[l], key, s[l - 1], t0, l);
		}
	}

	nmod_poly_clear(t0);
	nmod_poly_clear(beta);
	return result;
}

static int run(commit_t com[MSGS], nmod_poly_t m[MSGS], nmod_poly_t _m[MSGS],
		nmod_poly_t r[MSGS][WIDTH][2], commitkey_t *key, flint_rand_t rng) {
	int flag, result = 1;
	commit_t d[MSGS];
	nmod_poly_t t0, t1, rho, s[MSGS], u[MSGS][2];
	nmod_poly_t y[MSGS][WIDTH][2], _y[MSGS][WIDTH][2], t[MSGS][2], _t[MSGS][2];

	nmod_poly_init(t0, MODP);
	nmod_poly_init(t1, MODP);
	nmod_poly_init(rho, MODP);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(s[i], MODP);
		for (int k = 0; k < 2; k++) {
			nmod_poly_init(t[i][k], MODP);
			nmod_poly_init(_t[i][k], MODP);
			nmod_poly_init(u[i][k], MODP);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_init(y[i][j][k], MODP);
				nmod_poly_init(_y[i][j][k], MODP);
			}
		}
	}

	/* Verifier samples \rho that is different from the messages, and \beta. */
	do {
		flag = 1;
		commit_sample_rand(rho, rng);
		for (int i = 0; i < MSGS; i++) {
			if (nmod_poly_equal(rho, _m[i]) == 1) {
				flag = 0;
			}
		}
	} while (flag == 0);

	/* Verifier shifts the commitments by rho. */
	nmod_poly_rem(t0, rho, *commit_irred(0));
	nmod_poly_rem(t1, rho, *commit_irred(1));
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_sub(com[i].c2[0], com[i].c2[0], t0);
		nmod_poly_sub(com[i].c2[1], com[i].c2[1], t1);
	}

	shuffle_prover(y, _y, t, _t, u, d, s, com, m, _m, r, rho, key, rng);

	result = shuffle_verifier(y, _y, t, _t, u, d, s, com, _m, key);

	nmod_poly_clear(t0);
	nmod_poly_clear(t1);
	nmod_poly_clear(rho);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_clear(s[i]);
		for (int k = 0; k < 2; k++) {
			nmod_poly_clear(t[i][k]);
			nmod_poly_clear(_t[i][k]);
			nmod_poly_clear(u[i][k]);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_clear(y[i][j][k]);
				nmod_poly_clear(_y[i][j][k]);
			}
		}
	}

	return result;
}

static void test(flint_rand_t rand) {
	commitkey_t key;
	commit_t com[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS], r[MSGS][WIDTH][2];

	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < 2; k++) {
				nmod_poly_init(r[i][j][k], MODP);
			}
		}
	}

	/* Generate commitment key-> */
	commit_setup();
	commit_keygen(&key, rand);

	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short_crt(r[i][j]);
		}
		commit_sample_short(m[i]);
		commit_doit(&com[i], m[i], &key, r[i]);
	}

	/* Prover shuffles messages (only a circular shift for simplicity). */
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_set(_m[i], m[(i + 1) % MSGS]);
	}

	TEST_ONCE("shuffle proof is consistent") {
		TEST_ASSERT(run(com, m, _m, r, &key, rand) == 1, end);
	} TEST_END;

end:
	commit_finish();

	for (int i = 0; i < MSGS; i++) {
		commit_free(&com[i]);
		nmod_poly_clear(m[i]);
		nmod_poly_clear(_m[i]);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < 2; k++) {
				nmod_poly_clear(r[i][j][k]);
			}
		}
	}
	commit_keyfree(&key);
}

static void bench(flint_rand_t rand) {
	commitkey_t key;
	commit_t com[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS];
	nmod_poly_t alpha, beta, s[MSGS - 1];
	nmod_poly_t r[MSGS][WIDTH][2], y[WIDTH][2], _y[WIDTH][2];
	nmod_poly_t t[2], _t[2], u[2], v[2], _v[2];

	nmod_poly_init(alpha, MODP);
	nmod_poly_init(beta, MODP);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
		if (i != MSGS - 1) nmod_poly_init(s[i], MODP);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < 2; k++) {
				nmod_poly_init(r[i][j][k], MODP);
			}
		}
	}

	/* Generate commitment key-> */
	commit_setup();
	commit_keygen(&key, rand);

	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short_crt(r[i][j]);
		}
		commit_sample_short(m[i]);
		commit_doit(&com[i], m[i], &key, r[i]);
	}

	/* Prover shuffles messages (only a circular shift for simplicity). */
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_set(_m[i], m[(i + 1) % MSGS]);
	}

	BENCH_BEGIN("shuffle-proof (N messages)") {
		BENCH_ADD(run(com, m, _m, r, &key, rand));
	} BENCH_END;

	for (int i = 0; i < 2; i++) {
		nmod_poly_init(t[i], MODP);
		nmod_poly_init(_t[i], MODP);
		nmod_poly_init(u[i], MODP);
		nmod_poly_init(v[i], MODP);
		nmod_poly_init(_v[i], MODP);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_init(y[i][j], MODP);
			nmod_poly_init(_y[i][j], MODP);
		}
	}

	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short_crt(r[i][j]);
		}
	}
	commit_sample_rand(beta, rand);
	commit_sample_rand(alpha, rand);

	BENCH_BEGIN("linear-proof") {
		BENCH_ADD(lin_prover(y, _y, t, _t, u, com[0], com[1], &key, alpha, beta, r[0], r[0], 0));
	} BENCH_END;

	BENCH_BEGIN("verifier proof") {
		BENCH_ADD(lin_verifier(y, _y, t, _t, u, com[0], com[1], &key, alpha, beta, 0));
	} BENCH_END;

	commit_finish();

	nmod_poly_clear(alpha);
	nmod_poly_clear(beta);
	for (int i = 0; i < MSGS; i++) {
		commit_free(&com[i]);
		nmod_poly_clear(m[i]);
		nmod_poly_clear(_m[i]);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < 2; k++) {
				nmod_poly_clear(r[i][j][k]);
			}
		}
	}
	for (int i = 0; i < 2; i++) {
		nmod_poly_clear(t[i]);
		nmod_poly_clear(_t[i]);
		nmod_poly_clear(u[i]);
		nmod_poly_clear(v[i]);
		nmod_poly_clear(_v[i]);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_clear(y[i][j]);
			nmod_poly_clear(_y[i][j]);
		}
	}

	commit_keyfree(&key);
}

int main(int argc, char *argv[]) {
	flint_rand_t rand;

	flint_randinit(rand);

	printf("\n** Tests for lattice-based shuffle proof:\n\n");
	test(rand);

	printf("\n** Benchmarks for lattice-based shuffle proof:\n\n");
	bench(rand);

	flint_randclear(rand);
}
