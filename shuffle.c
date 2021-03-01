#include "param.h"
#include "commit.h"
#include "test.h"
#include "bench.h"
#include "assert.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#define MAXNORM 	44194

#define MSGS        2

static int run(flint_rand_t rand, commitkey_t *key) {
	int result = 1;
	commit_t com[2], d[2];
	nmod_poly_t tmp, m[2], _m[2], rho, beta, theta, s, _d[2];
	nmod_poly_t r[WIDTH][2], _r[WIDTH][2], y[WIDTH][2], _y[WIDTH][2];
	nmod_poly_t t[2], _t[2], u[2], v[2], _v[2], z[WIDTH], _z[WIDTH];

	nmod_poly_init(tmp, MODP);
	nmod_poly_init(rho, MODP);
	nmod_poly_init(beta, MODP);
	nmod_poly_init(theta, MODP);
	nmod_poly_init(s, MODP);
	for (int i = 0; i < 2; i++) {
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
		nmod_poly_init(t[i], MODP);
		nmod_poly_init(_t[i], MODP);
		nmod_poly_init(u[i], MODP);
		nmod_poly_init(v[i], MODP);
		nmod_poly_init(_v[i], MODP);
		nmod_poly_init(_d[i], MODP);
		nmod_poly_zero(t[i]);
		nmod_poly_zero(_t[i]);
		nmod_poly_zero(u[i]);
		nmod_poly_zero(v[i]);
		nmod_poly_zero(_v[i]);
	}
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(z[i], MODP);
		nmod_poly_init(_z[i], MODP);
		for (int j = 0; j < 2; j++) {
			nmod_poly_init(r[i][j], MODP);
			nmod_poly_init(_r[i][j], MODP);
			nmod_poly_init(y[i][j], MODP);
			nmod_poly_init(_y[i][j], MODP);
		}
	}

	for (int i = 0; i < WIDTH; i++) {
		commit_sample_short_crt(r[i]);
		commit_sample_short_crt(_r[i]);
	}

	commit_sample_short(m[0]);
	commit_doit(&com[0], m[0], key, r);
	commit_sample_short(m[1]);
	commit_doit(&com[1], m[1], key, r);

	/* Prover shuffles messages (only two in this case). */
	nmod_poly_set(_m[0], m[1]);
	nmod_poly_set(_m[1], m[0]);

	/* Verifier samples \rho that is different from the messages, and \beta. */
	do {
		uint32_t coeff;
		commit_sample_rand(rho, rand);
		for (int i = 0; i < DEGREE; i++) {
			coeff = nmod_poly_get_coeff_ui(rho, i) % MAXNORM;
			nmod_poly_set_coeff_ui(rho, i, coeff);
		}
	} while (nmod_poly_equal(rho, _m[0]) != 0 && nmod_poly_equal(rho, _m[1]));
	/* Check that norm is below float(sqrt(p)/(sqrt(2) - 1)). */
	assert(commit_norm_inf(rho) <= MAXNORM);
	commit_sample_short(beta);

	/* Prover shifts the messages by rho and shuffles them. */
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_sub(m[i], m[i], rho);
		nmod_poly_sub(_m[i], _m[i], rho);
	}

	/* Verifier shifts the commitment by rho. */
	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < 2; j++) {
			nmod_poly_rem(tmp, rho, *commit_irred(j));
			nmod_poly_sub(com[i].c2[j], com[i].c2[j], tmp);
		}
	}

	/* Prover samples theta_i and computes commitments D_i. */
	commit_sample_rand(theta, rand);
	nmod_poly_mulmod(tmp, theta, _m[0], *commit_poly());
	commit_doit(&d[0], tmp, key, r);
	nmod_poly_mulmod(tmp, theta, m[1], *commit_poly());
	commit_doit(&d[1], tmp, key, r);

	commit_sample_rand(beta, rand);
	nmod_poly_mulmod(s, theta, _m[0], *commit_poly());
	nmod_poly_mulmod(tmp, beta, m[0], *commit_poly());
	nmod_poly_sub(s, s, tmp);
	nmod_poly_invmod(tmp, _m[0], *commit_poly());
	nmod_poly_mulmod(s, s, tmp, *commit_poly());

	/* Now run \Prod_LIN instances, one for each commitment. */
	for (int l = 0; l < MSGS; l++) {
		for (int i = 0; i < WIDTH; i++) {
			commit_sample_gauss_crt(y[i]);
			commit_sample_gauss_crt(_y[i]);
		}
		for (int i = 0; i < 2; i++) {
			nmod_poly_zero(t[i]);
			nmod_poly_zero(_t[i]);
			nmod_poly_zero(v[i]);
			nmod_poly_zero(_v[i]);
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

		nmod_poly_zero(u[0]);
		nmod_poly_zero(u[1]);
		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < 2; j++) {
				nmod_poly_mulmod(tmp, key->b2[i][j], y[i][j], *commit_irred(j));
				if (l == 0) nmod_poly_mulmod(tmp, beta, tmp, *commit_irred(j));
				if (l == 1) nmod_poly_mulmod(tmp, s, tmp, *commit_irred(j));
				nmod_poly_add(u[j], u[j], tmp);
				nmod_poly_mulmod(tmp, key->b2[i][j], _y[i][j], *commit_irred(j));
				nmod_poly_sub(u[j], u[j], tmp);
			}
		}

		/* Verifier. */
		commit_sample_chall_crt(_d);

		/* Prover */
		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < 2; j++) {
				nmod_poly_mulmod(tmp, _d[j], r[i][j], *commit_irred(j));
				nmod_poly_add(y[i][j], y[i][j], tmp);
				nmod_poly_mulmod(tmp, _d[j], r[i][j], *commit_irred(j));
				nmod_poly_add(_y[i][j], _y[i][j], tmp);
			}
			pcrt_poly_rec(z[i], y[i]);
			pcrt_poly_rec(_z[i], _y[i]);
		}

		/* Verifier. */
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
		for (int i = 0; i < WIDTH; i++) {
			assert(commit_norm2_sqr(z[i]) <= (uint64_t)4 * DEGREE * SIGMA_E * SIGMA_E);
			assert(commit_norm2_sqr(_z[i]) <= (uint64_t)4 * DEGREE * SIGMA_E * SIGMA_E);
		}
		for (int j = 0; j < 2; j++) {
			nmod_poly_mulmod(tmp, _d[j], com[l].c1[j], *commit_irred(j));
			nmod_poly_add(t[j], t[j], tmp);
			nmod_poly_mulmod(tmp, _d[j], d[l].c1[j], *commit_irred(j));
			nmod_poly_add(_t[j], _t[j], tmp);
			result &= nmod_poly_equal(t[j], v[j]);
			result &= nmod_poly_equal(_t[j], _v[j]);
		}

		if (l == 0) {
			for (int j = 0; j < 2; j++) {
				nmod_poly_mulmod(t[j], beta, com[l].c2[j], *commit_irred(j));
				nmod_poly_mulmod(tmp, s, _m[l], *commit_poly());
				nmod_poly_add(t[j], t[j], tmp);
				nmod_poly_sub(t[j], t[j], d[l].c2[j]);
				nmod_poly_mulmod(t[j], t[j], _d[j], *commit_irred(j));
				nmod_poly_add(t[j], t[j], u[j]);
				nmod_poly_rem(t[j], t[j], *commit_irred(j));
			}
		} else {
			for (int j = 0; j < 2; j++) {
				nmod_poly_mulmod(t[j], s, com[l].c2[j], *commit_irred(j));
				nmod_poly_mulmod(tmp, beta, _m[l], *commit_poly());
				nmod_poly_add(t[j], t[j], tmp);
				nmod_poly_sub(t[j], t[j], d[l].c2[j]);
				nmod_poly_mulmod(t[j], t[j], _d[j], *commit_irred(j));
				nmod_poly_add(t[j], t[j], u[j]);
			}
		}

		nmod_poly_zero(u[0]);
		nmod_poly_zero(u[1]);
		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < 2; j++) {
				nmod_poly_mulmod(tmp, key->b2[i][j], y[i][j], *commit_irred(j));
				if (l == 0) nmod_poly_mulmod(tmp, beta, tmp, *commit_irred(j));
				if (l == 1) nmod_poly_mulmod(tmp, s, tmp, *commit_irred(j));
				nmod_poly_add(u[j], u[j], tmp);
				nmod_poly_mulmod(tmp, key->b2[i][j], _y[i][j], *commit_irred(j));
				nmod_poly_sub(u[j], u[j], tmp);
			}
		}
		for (int j = 0; j < 2; j++) {
			result &= nmod_poly_equal(t[j], u[j]);
		}
	}

	nmod_poly_clear(tmp);
	nmod_poly_clear(rho);
	nmod_poly_clear(beta);
	nmod_poly_clear(theta);
	nmod_poly_clear(s);
	for (int i = 0; i < 2; i++) {
		//commit_clear(&com[i]);
		nmod_poly_clear(m[i]);
		nmod_poly_clear(_m[i]);
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

	return result;
}

static void test(flint_rand_t rand) {
	commitkey_t key;

	commit_setup();
	/* Generate commitment key-> */
	commit_keygen(&key, rand);

	TEST_ONCE("shuffle proof is consistent") {
		TEST_ASSERT(run(rand, &key) == 1, end);
	} TEST_END;
end:
	commit_finish();

	commit_keyfree(&key);
}

static void bench(flint_rand_t rand) {
	commitkey_t key;

	commit_setup();
	/* Generate commitment key-> */
	commit_keygen(&key, rand);

	BENCH_BEGIN("shuffle-proof (2 messages)") {
		BENCH_ADD(run(rand, &key));
	} BENCH_END;

	commit_finish();

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
