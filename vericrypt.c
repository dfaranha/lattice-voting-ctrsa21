#include "param.h"
#include "commit.h"
#include "test.h"
#include "bench.h"
#include "encrypt.h"
#include "gaussian.h"
#include "fastrandombytes.h"
#include "sha.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

// Dimension of the message space in verifiable encryption.
#define VECTOR	3

/* Type that represents a verifiable encryption ciphertext. */
typedef struct _veritext_t {
	ciphertext_t cipher[VECTOR];
	fmpz_mod_poly_t c;
	fmpz_mod_poly_t r[VECTOR][DIM][2];
	fmpz_mod_poly_t e[VECTOR][DIM][2];
	fmpz_mod_poly_t e_[VECTOR][2];
	fmpz_mod_poly_t u[VECTOR];
} veritext_t;

static int vericrypt_test_norm(veritext_t *out) {
	int result;
	fmpz_t coeff, max, qdiv2, *q;
	fmpz_mod_poly_t t;
	fmpz_mod_ctx_t *ctx;

	fmpz_init(coeff);
	fmpz_init(max);
	fmpz_init(qdiv2);

	q = encrypt_large_modulus();
	ctx = encrypt_large_modulus_ctx();

	fmpz_mod_poly_init(t, *ctx);
	fmpz_set(qdiv2, *q);
	fmpz_divexact_ui(qdiv2, qdiv2, 2);
	fmpz_set_ui(max, 0);

	/* Compute norm_infty. */
	for (int i = 0; i < VECTOR; i++) {
		for (int j = 0; j < DIM; j++) {
			qcrt_poly_rec(t, out->r[i][j]);
			for (int k = 0; k < DEGREE; k++) {
				fmpz_mod_poly_get_coeff_fmpz(coeff, t, k, *ctx);
				if (fmpz_cmp(coeff, qdiv2) > 0)
					fmpz_sub(coeff, coeff, *q);
				fmpz_abs(coeff, coeff);
				if (fmpz_cmp(coeff, max) > 0)
					fmpz_set(max, coeff);
			}
			qcrt_poly_rec(t, out->e[i][j]);
			for (int k = 0; k < DEGREE; k++) {
				fmpz_mod_poly_get_coeff_fmpz(coeff, t, k, *ctx);
				if (fmpz_cmp(coeff, qdiv2) > 0)
					fmpz_sub(coeff, coeff, *q);
				fmpz_abs(coeff, coeff);
				if (fmpz_cmp(coeff, max) > 0)
					fmpz_set(max, coeff);
			}
		}
		qcrt_poly_rec(t, out->e_[i]);
		for (int k = 0; k < DEGREE; k++) {
			fmpz_mod_poly_get_coeff_fmpz(coeff, t, k, *ctx);
			if (fmpz_cmp(coeff, qdiv2) > 0)
				fmpz_sub(coeff, coeff, *q);
			fmpz_abs(coeff, coeff);
			if (fmpz_cmp(coeff, max) > 0)
				fmpz_set(coeff, max);
		}
	}
	for (int i = 0; i < VECTOR; i++) {
		for (int k = 0; k < DEGREE; k++) {
			fmpz_mod_poly_get_coeff_fmpz(coeff, out->u[i], k, *ctx);
			fmpz_abs(coeff, coeff);
			if (fmpz_cmp(coeff, max) > 0)
				fmpz_set(coeff, max);
		}
	}

	/* Actually compute sigma^2 to simplify comparison. */
	fmpz_set_ui(coeff, 6 * SIGMA_E);
	result = fmpz_cmp(max, coeff);

	fmpz_clear(coeff);
	fmpz_clear(max);
	fmpz_clear(qdiv2);
	return (result <= 0);
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void vericrypt_free(veritext_t *out) {
	for (int i = 0; i < VECTOR; i++) {
		fmpz_mod_poly_clear(out->u[i], *encrypt_modulus_ctx());
		for (int j = 0; j < DIM; j++) {
			fmpz_mod_poly_clear(out->e_[i][j], *encrypt_large_modulus_ctx());
			for (int k = 0; k < 2; k++) {
				fmpz_mod_poly_clear(out->r[i][j][k],
						*encrypt_large_modulus_ctx());
				fmpz_mod_poly_clear(out->e[i][j][k],
						*encrypt_large_modulus_ctx());
			}
		}
	}
}

void vericrypt_hash(uint8_t hash[SHA256HashSize], publickey_t *pk,
		fmpz_mod_poly_t t[VECTOR], fmpz_mod_poly_t u, veritext_t *out,
		ciphertext_t y[VECTOR], fmpz_mod_poly_t _u) {
	fmpz_poly_t s;
	fmpz_mod_ctx_t *ctx = encrypt_large_modulus_ctx();
	SHA256Context sha;
	char *str = NULL;

	fmpz_poly_init(s);
	SHA256Reset(&sha);

	/* Hash public key (A,t). */
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			for (int k = 0; k < 2; k++) {
				fmpz_mod_poly_get_fmpz_poly(s, pk->A[i][j][k], *ctx);
				str = fmpz_poly_get_str(s);
				SHA256Input(&sha, (const uint8_t *)str, strlen(str));
				free(str);
			}
			fmpz_mod_poly_get_fmpz_poly(s, pk->t[i][j], *ctx);
			str = fmpz_poly_get_str(s);
			SHA256Input(&sha, (const uint8_t *)str, strlen(str));
			free(str);
		}
	}
	/* Hash linear relation pair (T,u). */
	for (int i = 0; i < VECTOR; i++) {
		fmpz_mod_poly_get_fmpz_poly(s, t[i], *ctx);
		str = fmpz_poly_get_str(s);
		SHA256Input(&sha, (const uint8_t *)str, strlen(str));
		free(str);
	}
	fmpz_mod_poly_get_fmpz_poly(s, u, *ctx);
	str = fmpz_poly_get_str(s);
	SHA256Input(&sha, (const uint8_t *)str, strlen(str));
	free(str);

	/* Hash ciphertexts c = (v, w), y. */
	for (int i = 0; i < VECTOR; i++) {
		for (int j = 0; j < DIM; j++) {
			for (int k = 0; k < 2; k++) {
				fmpz_mod_poly_get_fmpz_poly(s, out->cipher[i].v[j][k], *ctx);
				str = fmpz_poly_get_str(s);
				SHA256Input(&sha, (const uint8_t *)str, strlen(str));
				free(str);
				fmpz_mod_poly_get_fmpz_poly(s, y[i].v[j][k], *ctx);
				str = fmpz_poly_get_str(s);
				SHA256Input(&sha, (const uint8_t *)str, strlen(str));
				free(str);
			}
		}
		for (int k = 0; k < 2; k++) {
			fmpz_mod_poly_get_fmpz_poly(s, out->cipher[i].w[k], *ctx);
			str = fmpz_poly_get_str(s);
			SHA256Input(&sha, (const uint8_t *)str, strlen(str));
			free(str);
			fmpz_mod_poly_get_fmpz_poly(s, y[i].w[k], *ctx);
			str = fmpz_poly_get_str(s);
			SHA256Input(&sha, (const uint8_t *)str, strlen(str));
			free(str);
		}
	}
	fmpz_mod_poly_get_fmpz_poly(s, _u, *encrypt_modulus_ctx());
	str = fmpz_poly_get_str(s);
	SHA256Input(&sha, (const uint8_t *)str, strlen(str));
	free(str);

	SHA256Result(&sha, hash);
	fmpz_poly_clear(s);
}

void vericrypt_sample_chall(fmpz_mod_poly_t f, uint8_t *seed, int len,
		fmpz_mod_ctx_t *ctx) {
	fmpz_t coeff;
	uint32_t buf;

	fmpz_init(coeff);
	fmpz_mod_poly_zero(f, *ctx);
	fastrandombytes_setseed(seed);
	for (int j = 0; j < NONZERO; j++) {
		fastrandombytes((unsigned char *)&buf, sizeof(buf));
		buf = buf % DEGREE;
		fmpz_mod_poly_get_coeff_fmpz(coeff, f, buf, *ctx);
		while (!fmpz_is_zero(coeff)) {
			fastrandombytes((unsigned char *)&buf, sizeof(buf));
			buf = buf % DEGREE;
			fmpz_mod_poly_get_coeff_fmpz(coeff, f, buf, *ctx);
		}
		fmpz_set_ui(coeff, 1);
		fmpz_mod_poly_set_coeff_fmpz(f, buf, coeff, *ctx);
	}
	fmpz_clear(coeff);
}

void vericrypt_sample_gauss(fmpz_mod_poly_t r, fmpz_mod_ctx_t *ctx) {
	fmpz_t coeff, t;
	fmpz_init(coeff);
	fmpz_init(t);
	fmpz_mod_poly_zero(r, *ctx);
	fmpz_mod_poly_fit_length(r, DEGREE, *ctx);
	for (int i = 0; i < DEGREE; i++) {
		int64_t sample = discrete_gaussian(0.0);
		fmpz_set_si(t, sample);
		fmpz_mod_poly_set_coeff_fmpz(r, i, t, *ctx);
	}
	fmpz_clear(coeff);
	fmpz_clear(t);
}

void vericrypt_sample_gauss_crt(fmpz_mod_poly_t r[2], fmpz_mod_ctx_t *ctx) {
	fmpz_mod_poly_t t;

	fmpz_mod_poly_init(t, *ctx);
	vericrypt_sample_gauss(t, ctx);
	fmpz_mod_poly_rem(r[0], t, *encrypt_irred(0), *ctx);
	fmpz_mod_poly_rem(r[1], t, *encrypt_irred(1), *ctx);

	fmpz_mod_poly_clear(t, *encrypt_large_modulus_ctx());
}

int vericrypt_doit(veritext_t *out, fmpz_mod_poly_t t[VECTOR],
		fmpz_mod_poly_t u, fmpz_mod_poly_t m[VECTOR], publickey_t *pk,
		flint_rand_t rand) {
	fmpz_mod_poly_t _u, tmp, c[DIM], y_mu[VECTOR];
	qcrt_poly_t y_r[VECTOR][DIM], y_e[VECTOR][DIM], y_e_[VECTOR];
	ciphertext_t y[VECTOR];
	uint8_t hash[SHA256HashSize];
	int result = 0;

	fmpz_mod_ctx_t *ctx_p = encrypt_modulus_ctx();
	fmpz_mod_ctx_t *ctx_q = encrypt_large_modulus_ctx();

	fmpz_mod_poly_init(tmp, *ctx_p);
	fmpz_mod_poly_init(_u, *ctx_p);
	fmpz_mod_poly_init(out->c, *ctx_p);
	fmpz_mod_poly_init(c[0], *ctx_p);
	fmpz_mod_poly_init(c[1], *ctx_p);
	for (int i = 0; i < VECTOR; i++) {
		for (int j = 0; j < DIM; j++) {
			for (int k = 0; k < 2; k++) {
				fmpz_mod_poly_init(y_r[i][j][k], *ctx_q);
				fmpz_mod_poly_init(y_e[i][j][k], *ctx_q);
				fmpz_mod_poly_init(out->r[i][j][k], *ctx_q);
				fmpz_mod_poly_init(out->e[i][j][k], *ctx_q);
			}
			fmpz_mod_poly_init(out->e_[i][j], *ctx_q);
			fmpz_mod_poly_init(y_e_[i][j], *ctx_q);
		}
		fmpz_mod_poly_init(out->u[i], *ctx_p);
		fmpz_mod_poly_init(y_mu[i], *ctx_p);
	}

	for (int i = 0; i < VECTOR; i++) {
		for (int j = 0; j < DIM; j++) {
			encrypt_sample_short_crt(out->r[i][j], *ctx_q);
			encrypt_sample_short_crt(out->e[i][j], *ctx_q);
		}
		encrypt_sample_short_crt(out->e_[i], *ctx_q);
	}

	for (int i = 0; i < VECTOR; i++) {
		encrypt_make(&out->cipher[i], out->r[i], out->e[i], out->e_[i], m[i],
				pk);
	}

	while (result == 0) {
		for (int i = 0; i < VECTOR; i++) {
			for (int j = 0; j < DIM; j++) {
				vericrypt_sample_gauss_crt(y_r[i][j], ctx_q);
				vericrypt_sample_gauss_crt(y_e[i][j], ctx_q);
			}
			vericrypt_sample_gauss_crt(y_e_[i], ctx_q);
			vericrypt_sample_gauss(y_mu[i], ctx_p);
		}

		for (int i = 0; i < VECTOR; i++) {
			encrypt_make(&y[i], y_r[i], y_e[i], y_e_[i], y_mu[i], pk);
		}

		fmpz_mod_poly_init(_u, *ctx_p);
		for (int i = 0; i < VECTOR; i++) {
			fmpz_mod_poly_mulmod(tmp, t[i], y_mu[i], *encrypt_poly(), *ctx_p);
			fmpz_mod_poly_add(_u, _u, tmp, *ctx_p);
		}

		/* Hash and convert result to challenge space. */
		vericrypt_hash(hash, pk, t, u, out, y, _u);
		vericrypt_sample_chall(out->c, hash, SHA256HashSize, ctx_p);
		fmpz_mod_poly_rem(c[0], out->c, *encrypt_irred(0), *ctx_q);
		fmpz_mod_poly_rem(c[1], out->c, *encrypt_irred(1), *ctx_q);

		// Compute z = [r e e' mu]^T c + y.
		for (int i = 0; i < VECTOR; i++) {
			for (int j = 0; j < DIM; j++) {
				for (int k = 0; k < 2; k++) {
					fmpz_mod_poly_mulmod(out->r[i][j][k], out->r[i][j][k], c[k],
							*encrypt_irred(k), *ctx_q);
					fmpz_mod_poly_add(out->r[i][j][k], out->r[i][j][k],
							y_r[i][j][k], *ctx_q);
					fmpz_mod_poly_mulmod(out->e[i][j][k], out->e[i][j][k], c[k],
							*encrypt_irred(k), *ctx_q);
					fmpz_mod_poly_add(out->e[i][j][k], out->e[i][j][k],
							y_e[i][j][k], *ctx_q);
				}
			}
		}
		for (int i = 0; i < VECTOR; i++) {
			for (int k = 0; k < 2; k++) {
				fmpz_mod_poly_mulmod(out->e_[i][k], out->e_[i][k], c[k],
						*encrypt_irred(k), *ctx_q);
				fmpz_mod_poly_add(out->e_[i][k], out->e_[i][k], y_e_[i][k],
						*ctx_q);
			}
			fmpz_mod_poly_mulmod(out->u[i], m[i], out->c, *encrypt_poly(),
					*ctx_p);
			fmpz_mod_poly_add(out->u[i], out->u[i], y_mu[i], *ctx_p);
		}

		result = vericrypt_test_norm(out);
	}

	fmpz_mod_poly_clear(tmp, *ctx_p);
	fmpz_mod_poly_clear(c[0], *ctx_p);
	fmpz_mod_poly_clear(c[1], *ctx_p);
	for (int i = 0; i < VECTOR; i++) {
		for (int j = 0; j < DIM; j++) {
			for (int k = 0; k < 2; k++) {
				fmpz_mod_poly_clear(y_r[i][j][k], *ctx_q);
				fmpz_mod_poly_clear(y_e[i][j][k], *ctx_q);
			}
			fmpz_mod_poly_clear(y_e_[i][j], *ctx_q);
		}
		fmpz_mod_poly_clear(y_mu[i], *ctx_p);
	}
	return result;
}

int vericrypt_verify(veritext_t *in, fmpz_mod_poly_t t[VECTOR],
		fmpz_mod_poly_t u, publickey_t *pk) {
	fmpz_mod_poly_t _u, tp, tq, c;
	qcrt_poly_t _c;
	ciphertext_t y[VECTOR];
	int result = 0;
	uint8_t hash[SHA256HashSize];

	fmpz_mod_ctx_t *ctx_p = encrypt_modulus_ctx();
	fmpz_mod_ctx_t *ctx_q = encrypt_large_modulus_ctx();

	fmpz_mod_poly_init(c, *ctx_p);
	fmpz_mod_poly_init(tq, *ctx_q);
	fmpz_mod_poly_init(tp, *ctx_p);
	fmpz_mod_poly_init(_u, *ctx_p);
	for (int i = 0; i < DIM; i++) {
		fmpz_mod_poly_init(_c[i], *ctx_q);
	}

	if (vericrypt_test_norm(in)) {
		for (int i = 0; i < VECTOR; i++) {
			encrypt_make(&y[i], in->r[i], in->e[i], in->e_[i], in->u[i], pk);
		}
		fmpz_mod_poly_rem(_c[0], in->c, *encrypt_irred(0), *ctx_q);
		fmpz_mod_poly_rem(_c[1], in->c, *encrypt_irred(1), *ctx_q);

		fmpz_mod_poly_zero(_u, *ctx_p);
		for (int i = 0; i < VECTOR; i++) {
			fmpz_mod_poly_mulmod(tp, t[i], in->u[i], *encrypt_poly(), *ctx_p);
			fmpz_mod_poly_add(_u, _u, tp, *ctx_p);
		}
		fmpz_mod_poly_mulmod(tp, in->c, u, *encrypt_poly(), *ctx_p);
		fmpz_mod_poly_sub(_u, _u, tp, *ctx_p);

		for (int i = 0; i < VECTOR; i++) {
			for (int j = 0; j < DIM; j++) {
				for (int k = 0; k < 2; k++) {
					fmpz_mod_poly_mulmod(tq, _c[k], in->cipher[i].v[j][k],
							*encrypt_irred(k), *ctx_q);
					fmpz_mod_poly_sub(y[i].v[j][k], y[i].v[j][k], tq, *ctx_q);
				}
			}
			for (int k = 0; k < 2; k++) {
				fmpz_mod_poly_mulmod(tq, _c[k], in->cipher[i].w[k],
						*encrypt_irred(k), *ctx_q);
				fmpz_mod_poly_sub(y[i].w[k], y[i].w[k], tq, *ctx_q);
			}
		}

		vericrypt_hash(hash, pk, t, u, in, y, _u);
		vericrypt_sample_chall(c, hash, SHA256HashSize, ctx_p);

		result = fmpz_mod_poly_equal(in->c, c, *ctx_p);
	}

	fmpz_mod_poly_clear(c, *ctx_p);
	fmpz_mod_poly_clear(tp, *ctx_p);
	fmpz_mod_poly_clear(tq, *ctx_q);
	for (int i = 0; i < DIM; i++) {
		fmpz_mod_poly_clear(_c[i], *ctx_q);
	}
	return (result == 1);
}

int vericrypt_undo(fmpz_mod_poly_t m[VECTOR], veritext_t *in,
		fmpz_mod_poly_t t[VECTOR], fmpz_mod_poly_t u, publickey_t *pk,
		privatekey_t *sk) {
	fmpz_mod_poly_t _c;

	if (vericrypt_verify(in, t, u, pk) == 0) {
		return 0;
	}

	fmpz_mod_poly_init(_c, *encrypt_modulus_ctx());
	encrypt_sample_short(_c, *encrypt_modulus_ctx());
	fmpz_mod_poly_sub(_c, in->c, _c, *encrypt_modulus_ctx());

	for (int i = 0; i < VECTOR; i++) {
		encrypt_undo(m[i], _c, &in->cipher[i], sk);
	}

	fmpz_mod_poly_clear(_c, *encrypt_modulus_ctx());
	return 1;
}

#ifdef MAIN
// Tests and benchmarks below.
static void test(flint_rand_t rand) {
	publickey_t pk;
	privatekey_t sk;
	veritext_t c;
	fmpz_mod_poly_t tmp, t[VECTOR], u;
	fmpz_mod_poly_t m[VECTOR], _m[VECTOR];

	fmpz_mod_poly_init(tmp, *encrypt_modulus_ctx());
	for (int i = 0; i < VECTOR; i++) {
		fmpz_mod_poly_init(m[i], *encrypt_modulus_ctx());
		fmpz_mod_poly_init(_m[i], *encrypt_modulus_ctx());
		encrypt_sample_short(m[i], *encrypt_modulus_ctx());
	}
	for (int i = 0; i < VECTOR; i++) {
		fmpz_mod_poly_init(t[i], *encrypt_modulus_ctx());
	}
	fmpz_mod_poly_init(u, *encrypt_modulus_ctx());

	encrypt_keygen(&pk, &sk, rand);

	TEST_BEGIN("verifiable encryption is consistent") {
		for (int i = 0; i < VECTOR; i++) {
			encrypt_sample_short(m[i], *encrypt_modulus_ctx());
		}

		fmpz_mod_poly_zero(u, *encrypt_modulus_ctx());
		for (int i = 0; i < VECTOR; i++) {
			fmpz_mod_poly_randtest(t[i], rand, DEGREE, *encrypt_modulus_ctx());
			fmpz_mod_poly_mulmod(tmp, t[i], m[i], *encrypt_poly(),
					*encrypt_modulus_ctx());
			fmpz_mod_poly_add(u, u, tmp, *encrypt_modulus_ctx());
		}

		TEST_ASSERT(vericrypt_doit(&c, t, u, m, &pk, rand) == 1, end);
		TEST_ASSERT(vericrypt_verify(&c, t, u, &pk) == 1, end);
		TEST_ASSERT(vericrypt_undo(_m, &c, t, u, &pk, &sk) == 1, end);

		for (int i = 0; i < VECTOR; i++) {
			//TEST_ASSERT(fmpz_mod_poly_equal(m[i], _m[i], *encrypt_modulus_ctx()) == 1, end);
		}
	} TEST_END;

  end:
	encrypt_keyfree(&pk, &sk);

	fmpz_mod_poly_clear(tmp, *encrypt_modulus_ctx());
	for (int i = 0; i < VECTOR; i++) {
		fmpz_mod_poly_clear(m[i], *encrypt_modulus_ctx());
		fmpz_mod_poly_clear(_m[i], *encrypt_modulus_ctx());
	}
	for (int i = 0; i < VECTOR; i++) {
		fmpz_mod_poly_clear(t[i], *encrypt_modulus_ctx());
	}
	fmpz_mod_poly_clear(u, *encrypt_modulus_ctx());
}

static void bench(flint_rand_t rand) {
	publickey_t pk;
	privatekey_t sk;
	veritext_t c;
	fmpz_mod_poly_t tmp, t[VECTOR], u;
	fmpz_mod_poly_t m[VECTOR], _m[VECTOR];

	fmpz_mod_poly_init(tmp, *encrypt_modulus_ctx());
	for (int i = 0; i < VECTOR; i++) {
		fmpz_mod_poly_init(m[i], *encrypt_modulus_ctx());
		fmpz_mod_poly_init(_m[i], *encrypt_modulus_ctx());
		encrypt_sample_short(m[i], *encrypt_modulus_ctx());
	}
	for (int i = 0; i < VECTOR; i++) {
		fmpz_mod_poly_init(t[i], *encrypt_modulus_ctx());
	}
	fmpz_mod_poly_init(u, *encrypt_modulus_ctx());

	encrypt_keygen(&pk, &sk, rand);

	fmpz_mod_poly_zero(u, *encrypt_modulus_ctx());
	for (int i = 0; i < VECTOR; i++) {
		fmpz_mod_poly_randtest(t[i], rand, DEGREE, *encrypt_modulus_ctx());
		fmpz_mod_poly_mulmod(tmp, t[i], m[i], *encrypt_poly(),
				*encrypt_modulus_ctx());
		fmpz_mod_poly_add(u, u, tmp, *encrypt_modulus_ctx());
	}

	BENCH_BEGIN("vericrypt_doit") {
		BENCH_ADD(vericrypt_doit(&c, t, u, m, &pk, rand));
	} BENCH_END;

	BENCH_BEGIN("vericrypt_verify") {
		BENCH_ADD(vericrypt_verify(&c, t, u, &pk));
	} BENCH_END;

	BENCH_BEGIN("vericrypt_sample_gauss") {
		BENCH_ADD(vericrypt_sample_gauss(tmp, encrypt_large_modulus_ctx()));
	} BENCH_END;

	BENCH_BEGIN("vericrypt_undo") {
		BENCH_ADD(vericrypt_undo(_m, &c, t, u, &pk, &sk));
	} BENCH_END;

	encrypt_keyfree(&pk, &sk);

	fmpz_mod_poly_clear(tmp, *encrypt_modulus_ctx());
	for (int i = 0; i < VECTOR; i++) {
		fmpz_mod_poly_clear(m[i], *encrypt_modulus_ctx());
		fmpz_mod_poly_clear(_m[i], *encrypt_modulus_ctx());
	}
	for (int i = 0; i < VECTOR; i++) {
		fmpz_mod_poly_clear(t[i], *encrypt_modulus_ctx());
	}
	fmpz_mod_poly_clear(u, *encrypt_modulus_ctx());
}

int main(int argc, char *argv[]) {
	flint_rand_t rand;

	flint_randinit(rand);
	encrypt_setup();

	printf("\n** Tests for lattice-based verifiable encryption:\n\n");
	test(rand);

	printf("\n** Benchmarks for lattice-based verifiable encryption:\n\n");
	bench(rand);

	encrypt_finish();
}

#endif
