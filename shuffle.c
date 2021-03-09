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

#define MSGS        2

void shuffle_hash(nmod_poly_t d[2], commitkey_t *key, commit_t x, commit_t y,
		nmod_poly_t alpha, nmod_poly_t beta, nmod_poly_t u[2],
		nmod_poly_t t[2], nmod_poly_t _t[2]) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	nmod_poly_t f;
	uint32_t buf;

	nmod_poly_init(f, MODP);
	SHA256Reset(&sha);

	/* Hash public key. */
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < 2; k++) {
				SHA256Input(&sha, (const uint8_t*)key->B1[i][j][k]->coeffs, key->B1[i][j][k]->alloc * sizeof(uint64_t));
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
	nmod_poly_zero(f);

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
	nmod_poly_sub(f, d[0], d[1]);
	nmod_poly_rem(d[0], f, *commit_irred(0));
	nmod_poly_rem(d[1], f, *commit_irred(1));
	nmod_poly_clear(f);
}

static void prover_lin(nmod_poly_t y[WIDTH][2], nmod_poly_t _y[WIDTH][2],
		nmod_poly_t t[2], nmod_poly_t _t[2], nmod_poly_t u[2],
		commit_t x, commit_t _x, commitkey_t *key,
		nmod_poly_t alpha, nmod_poly_t beta, nmod_poly_t r[WIDTH][2], int l) {
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
	shuffle_hash(d, key, x, _x, alpha, beta, u, t, _t);

	/* Prover */
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_mulmod(tmp, d[j], r[i][j], *commit_irred(j));
			nmod_poly_add(y[i][j], y[i][j], tmp);
			nmod_poly_mulmod(tmp, d[j], r[i][j], *commit_irred(j));
			nmod_poly_add(_y[i][j], _y[i][j], tmp);
		}
	}

	nmod_poly_clear(tmp);
	for (int i = 0; i < 2; i++) {
		nmod_poly_clear(d[i]);
	}
}

static int verifier_lin(commit_t com, commit_t x,
		nmod_poly_t _m[MSGS], nmod_poly_t y[WIDTH][2], nmod_poly_t _y[WIDTH][2],
		nmod_poly_t t[2], nmod_poly_t _t[2], nmod_poly_t u[2], commitkey_t *key,
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
	shuffle_hash(_d, key, com, x, alpha, beta, u, t, _t);

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

static int run(commit_t com[MSGS], nmod_poly_t m[MSGS], nmod_poly_t _m[MSGS],
		nmod_poly_t r[WIDTH][2], commitkey_t *key, flint_rand_t rng) {
	int flag, result = 1;
	commit_t d[MSGS];
	nmod_poly_t t0, t1, rho, beta, theta[MSGS], s[MSGS], _d[MSGS];
	nmod_poly_t y[WIDTH][2], _y[WIDTH][2], t[2], _t[2], u[2];

	nmod_poly_init(t0, MODP);
	nmod_poly_init(t1, MODP);
	nmod_poly_init(rho, MODP);
	nmod_poly_init(beta, MODP);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(theta[i], MODP);
		nmod_poly_init(s[i], MODP);
	}
	for (int i = 0; i < 2; i++) {
		nmod_poly_init(t[i], MODP);
		nmod_poly_init(_t[i], MODP);
		nmod_poly_init(u[i], MODP);
		nmod_poly_init(_d[i], MODP);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_init(y[i][j], MODP);
			nmod_poly_init(_y[i][j], MODP);
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

	/* Prover shifts the messages by rho and shuffles them. */
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_sub(m[i], m[i], rho);
		nmod_poly_sub(_m[i], _m[i], rho);
	}

	/* Verifier shifts the commitment by rho. */
	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_rem(t0, rho, *commit_irred(j));
			nmod_poly_sub(com[i].c2[j], com[i].c2[j], t0);
		}
	}

	/* Prover samples theta_i and computes commitments D_i. */
	commit_sample_rand(theta[0], rng);
	nmod_poly_mulmod(t0, theta[0], _m[0], *commit_poly());
	commit_doit(&d[0], t0, key, r);
	for (int i = 1; i < MSGS - 1; i++) {
		commit_sample_rand(theta[i], rng);
		nmod_poly_mulmod(t0, theta[i - 1], m[i], *commit_poly());
		nmod_poly_mulmod(t1, theta[i], _m[i], *commit_poly());
		nmod_poly_add(t0, t0, t1);
		commit_doit(&d[i], t0, key, r);
	}
	nmod_poly_mulmod(t0, theta[MSGS - 2], m[MSGS - 1], *commit_poly());
	commit_doit(&d[MSGS - 1], t0, key, r);

	commit_sample_rand(beta, rng);
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
	if (MSGS > 2) {
		nmod_poly_mulmod(s[MSGS - 1], theta[MSGS - 2], m[MSGS - 1], *commit_poly());
		nmod_poly_mulmod(t0, beta, _m[MSGS - 1], *commit_poly());
		if (MSGS & 1) {
			nmod_poly_sub(s[MSGS - 1], s[0], t0);
		} else {
			nmod_poly_add(s[MSGS - 1], s[0], t0);
		}
		nmod_poly_invmod(t0, m[MSGS - 1], *commit_poly());
		nmod_poly_mulmod(s[MSGS - 1], s[MSGS - 1], t0, *commit_poly());
	}

	/* Now run \Prod_LIN instances, one for each commitment. */
	for (int l = 0; l < MSGS; l++) {
		if (l < MSGS - 1) {
			nmod_poly_mulmod(t0, s[l], _m[l], *commit_poly());
		} else {
			nmod_poly_mulmod(t0, beta, _m[l], *commit_poly());
		}

		if (l == 0) {
			prover_lin(y, _y, t, _t, u, com[0], d[0], key, beta, t0, r, l);
			result &= verifier_lin(com[0], d[0], _m, y, _y, t, _t, u, key, beta, t0, l);
		} else {
			prover_lin(y, _y, t, _t, u, com[l], d[l], key, s[l - 1], t0, r, l);
			result &= verifier_lin(com[l], d[l], _m, y, _y, t, _t, u, key, s[l - 1], t0, l);
		}
	}

	nmod_poly_clear(t0);
	nmod_poly_clear(t1);
	nmod_poly_clear(rho);
	nmod_poly_clear(beta);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_clear(theta[i]);
		nmod_poly_clear(s[i]);
	}
	for (int i = 0; i < 2; i++) {
		nmod_poly_clear(t[i]);
		nmod_poly_clear(_t[i]);
		nmod_poly_clear(u[i]);
		nmod_poly_clear(_d[i]);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_clear(y[i][j]);
			nmod_poly_clear(_y[i][j]);
		}
	}

	return result;
}

static void test(flint_rand_t rand) {
	commitkey_t key;
	commit_t com[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS], r[WIDTH][2];

	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_init(r[i][j], MODP);
		}
	}

	/* Generate commitment key-> */
	commit_setup();
	commit_keygen(&key, rand);

	for (int i = 0; i < WIDTH; i++) {
		commit_sample_short_crt(r[i]);
	}
	for (int i = 0; i < MSGS; i++) {
		commit_sample_short(m[i]);
		commit_doit(&com[i], m[i], &key, r);
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
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_clear(r[i][j]);
		}
	}
	commit_keyfree(&key);
}

static void bench(flint_rand_t rand) {
	commitkey_t key;
	commit_t com[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS];
	nmod_poly_t alpha, beta, s[MSGS - 1], _d[MSGS];
	nmod_poly_t r[WIDTH][2], y[WIDTH][2], _y[WIDTH][2];
	nmod_poly_t t[2], _t[2], u[2], v[2], _v[2], z[WIDTH], _z[WIDTH];

	nmod_poly_init(alpha, MODP);
	nmod_poly_init(beta, MODP);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
		if (i != MSGS - 1) nmod_poly_init(s[i], MODP);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_init(r[i][j], MODP);
		}
	}

	/* Generate commitment key-> */
	commit_setup();
	commit_keygen(&key, rand);

	for (int i = 0; i < WIDTH; i++) {
		commit_sample_short_crt(r[i]);
	}
	for (int i = 0; i < MSGS; i++) {
		commit_sample_short(m[i]);
		commit_doit(&com[i], m[i], &key, r);
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
		nmod_poly_init(_d[i], MODP);
	}
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(z[i], MODP);
		nmod_poly_init(_z[i], MODP);
		for (int j = 0; j < 2; j++) {
			nmod_poly_init(r[i][j], MODP);
			nmod_poly_init(y[i][j], MODP);
			nmod_poly_init(_y[i][j], MODP);
		}
	}

	for (int i = 0; i < WIDTH; i++) {
		commit_sample_short_crt(r[i]);
	}
	commit_sample_rand(beta, rand);
	commit_sample_rand(alpha, rand);

	BENCH_BEGIN("linear-proof") {
		BENCH_ADD(prover_lin(y, _y, t, _t, u, com[0], com[1], &key, alpha, beta, r, 0));
	} BENCH_END;

	BENCH_BEGIN("verifier proof") {
		BENCH_ADD(verifier_lin(com[0], com[1], _m, y, _y, t, _t, u, &key, alpha, beta, 0));
	} BENCH_END;

	commit_finish();

	nmod_poly_clear(alpha);
	nmod_poly_clear(beta);
	for (int i = 0; i < MSGS; i++) {
		commit_free(&com[i]);
		nmod_poly_clear(m[i]);
		nmod_poly_clear(_m[i]);
	}
	for (int i = 0; i < 2; i++) {
		nmod_poly_clear(t[i]);
		nmod_poly_clear(_t[i]);
		nmod_poly_clear(u[i]);
		nmod_poly_clear(v[i]);
		nmod_poly_clear(_v[i]);
		nmod_poly_clear(_d[i]);

	}
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_clear(z[i]);
		nmod_poly_clear(_z[i]);
		for (int j = 0; j < 2; j++) {
			nmod_poly_clear(r[i][j]);
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
}
