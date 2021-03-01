/**
 * @file
 *
 * Implementation of the lattice-based encryption scheme.
 *
 * @ingroup commit
 */

#include "param.h"
#include "test.h"
#include "bench.h"
#include "encrypt.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/* The large modulus for the encryption scheme. */
#define Q		"72057594037928893"
#define Q0		"29973109198516688"
#define Q1		"42084484839412205"
#define Q2		"36028797018964446"

/* Prime modulus for defining commitment ring. */
static fmpz_t p;

/* Prime modulus for defining encryption ring. */
static fmpz_t q;

/** Context for arithmetic modulo q. */
static fmpz_mod_ctx_t ctx_q;

/** Context for arithmetic modulo p. */
static fmpz_mod_ctx_t ctx_p;

/* Polynomial defining the cyclotomic ring. */
static fmpz_mod_poly_t large_poly, poly;

/* Pairs of irreducible polynomials for CRT representation. */
static qcrt_poly_t irred;

/* Inverses of the irreducible polynomials for CRT reconstruction. */
static qcrt_poly_t inv;

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

/* Recover polynomial from CRT representation. */
void qcrt_poly_rec(fmpz_mod_poly_t c, qcrt_poly_t a) {
	fmpz_mod_poly_t t;

	fmpz_mod_poly_init(t, ctx_q);

	fmpz_mod_poly_sub(t, a[0], a[1], ctx_q);
	fmpz_mod_poly_mul(t, t, inv[1], ctx_q);
	fmpz_mod_poly_mul(c, t, irred[1], ctx_q);
	fmpz_mod_poly_add(c, c, a[1], ctx_q);
	fmpz_mod_poly_mul(t, irred[0], irred[1], ctx_q);
	fmpz_mod_poly_rem(c, c, t, ctx_q);

	fmpz_mod_poly_clear(t, ctx_q);
}

// Sample short element.
void encrypt_sample_short(fmpz_mod_poly_t r, fmpz_mod_ctx_t ctx) {
	uint64_t buf;
	fmpz_t coeff;

	fmpz_init(coeff);
	fmpz_mod_poly_zero(r, ctx);
	fmpz_mod_poly_fit_length(r, DEGREE, ctx);
	for (int i = 0; i < DEGREE; i += 32) {
		getrandom(&buf, sizeof(buf), 0);
		for (int j = 0; j < 64; j += 2) {
			fmpz_set_ui(coeff, ((buf >> (j + 1)) & 1));
			if ((buf >> j) & 1) {
				fmpz_neg(coeff, coeff);
			}
			fmpz_mod_poly_set_coeff_fmpz(r, (i+j/2) % DEGREE, coeff, ctx);
		}
	}

	fmpz_clear(coeff);
}

// Sample short element in CRT representation.
void encrypt_sample_short_crt(fmpz_mod_poly_t r[2], fmpz_mod_ctx_t ctx) {
	fmpz_mod_poly_t t;

	fmpz_mod_poly_init(t, ctx);
	encrypt_sample_short(t, ctx);
	fmpz_mod_poly_rem(r[0], t, irred[0], ctx);
	fmpz_mod_poly_rem(r[1], t, irred[1], ctx);

	fmpz_mod_poly_clear(t, ctx);
}

// Initialize encryption scheme.
void encrypt_setup() {
	fmpz_t q0, q1;

	fmpz_init(p);
	fmpz_init(q);
	fmpz_init(q0);
	fmpz_init(q1);

	fmpz_set_ui(p, MODP);
	fmpz_set_str(q, Q, 10);
	fmpz_set_str(q0, Q0, 10);
	fmpz_set_str(q1, Q1, 10);

	fmpz_mod_ctx_init(ctx_p, p);
	fmpz_mod_ctx_init(ctx_q, q);

	fmpz_mod_poly_init(poly, ctx_p);
	fmpz_mod_poly_init(large_poly, ctx_q);
	for (int i = 0; i < 2; i++) {
		fmpz_mod_poly_init(irred[i], ctx_q);
		fmpz_mod_poly_init(inv[i], ctx_q);
	}

	// Initialize cyclotomic polynomial (x^N + 1) over F_p
	fmpz_mod_poly_set_coeff_ui(poly, DEGREE, 1, ctx_p);
	fmpz_mod_poly_set_coeff_ui(poly, 0, 1, ctx_p);

	// Initialize cyclotomic polynomial (x^N + 1) over F_q
	fmpz_mod_poly_set_coeff_ui(large_poly, DEGREE, 1, ctx_q);
	fmpz_mod_poly_set_coeff_ui(large_poly, 0, 1, ctx_q);

	// Initialize each factor as well.
	fmpz_mod_poly_set_coeff_ui(irred[0], DEGCRT, 1, ctx_q);
	fmpz_mod_poly_set_coeff_fmpz(irred[0], 0, q0, ctx_q);
	fmpz_mod_poly_set_coeff_ui(irred[1], DEGCRT, 1, ctx_q);
	fmpz_mod_poly_set_coeff_fmpz(irred[1], 0, q1, ctx_q);

	fmpz_mod_poly_invmod(inv[0], irred[0], irred[1], ctx_q);
	fmpz_mod_poly_invmod(inv[1], irred[1], irred[0], ctx_q);

	fmpz_clear(q0);
	fmpz_clear(q1);
}

// Return small modulus p.
fmpz_t *encrypt_modulus() {
	return &p;
}

// Return large modulus q.
fmpz_t *encrypt_large_modulus() {
	return &q;
}

// Return small modulus p.
fmpz_mod_ctx_t *encrypt_modulus_ctx() {
	return &ctx_p;
}

// Return large modulus q.
fmpz_mod_ctx_t *encrypt_large_modulus_ctx() {
	return &ctx_q;
}

// Return cyclotomic polynomial.
fmpz_mod_poly_t *encrypt_large_poly() {
	return &large_poly;
}

// Return cyclotomic polynomial.
fmpz_mod_poly_t *encrypt_poly() {
	return &poly;
}

// Return irreducible polynomials for CRT representation.
fmpz_mod_poly_t *encrypt_irred(int i) {
	return &irred[i];
}

// Finalize encryption scheme.
void encrypt_finish() {
	fmpz_mod_poly_clear(poly, ctx_p);
	fmpz_mod_poly_clear(large_poly, ctx_q);
	for (int i = 0; i < 2; i++) {
		fmpz_mod_poly_clear(irred[i], ctx_q);
		fmpz_mod_poly_clear(inv[i], ctx_q);
	}
	fmpz_mod_ctx_clear(ctx_p);
	fmpz_mod_ctx_clear(ctx_q);
	fmpz_clear(p);
	fmpz_clear(q);
}

// Generate a key pair.
void encrypt_keygen(publickey_t *pk, privatekey_t *sk, flint_rand_t rand) {
	fmpz_mod_poly_t t;

	fmpz_mod_poly_init(t, ctx_q);
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < 2; j++) {
			fmpz_mod_poly_init(sk->s1[i][j], ctx_q);
			fmpz_mod_poly_init(sk->s2[i][j], ctx_q);
		}
		encrypt_sample_short_crt(sk->s1[i], ctx_q);
		encrypt_sample_short_crt(sk->s2[i], ctx_q);
	}
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			fmpz_mod_poly_init(pk->t[i][j], ctx_q);
			fmpz_mod_poly_zero(pk->t[i][j], ctx_q);
			for (int k = 0; k < 2; k++) {
				fmpz_mod_poly_init(pk->A[i][j][k], ctx_q);
				fmpz_mod_poly_randtest(pk->A[i][j][k], rand, DEGCRT, ctx_q);
			}
		}
		for (int k = 0; k < 2; k++) {
			fmpz_mod_poly_add(pk->t[i][k], pk->t[i][k], sk->s2[i][k], ctx_q);
		}
	}

	// Compute (A, t = As_1 + s_2).
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			for (int k = 0; k < 2; k++) {
				fmpz_mod_poly_mulmod(t, pk->A[i][j][k], sk->s1[j][k], irred[k], ctx_q);
				fmpz_mod_poly_add(pk->t[i][k], pk->t[i][k], t, ctx_q);
			}
		}
	}
	fmpz_mod_poly_clear(t, ctx_q);
}

// Free key pair.
void encrypt_keyfree(publickey_t *pk, privatekey_t *sk) {
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			fmpz_mod_poly_clear(pk->t[i][j], ctx_q);
			for (int k = 0; k < 2; k++) {
				fmpz_mod_poly_clear(pk->A[i][j][k], ctx_q);
			}
		}
		for (int k = 0; k < 2; k++) {
			fmpz_mod_poly_clear(sk->s1[i][k], ctx_q);
			fmpz_mod_poly_clear(sk->s2[i][k], ctx_q);
		}
	}
}

// Internal encryption function.
void encrypt_make(ciphertext_t *c, qcrt_poly_t r[DIM], qcrt_poly_t e[DIM],
	qcrt_poly_t e_, fmpz_mod_poly_t m, publickey_t *pk) {
	fmpz_poly_t s;
	fmpz_mod_poly_t _m, t;
	fmpz_t coeff, p2;

	fmpz_init(coeff);
	fmpz_init(p2);

	fmpz_poly_init(s);
	fmpz_mod_poly_init(_m, ctx_q);
	for (int i = 0; i < DIM; i++) {
		fmpz_mod_poly_init(c->w[i], ctx_q);
		for (int j = 0; j < 2; j++) {
			fmpz_mod_poly_init(c->v[i][j], ctx_q);
			fmpz_mod_poly_zero(c->v[i][j], ctx_q);
		}
	}

	fmpz_mod_poly_init(t, ctx_q);
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			for (int k = 0; k < 2; k++) {
				fmpz_mod_poly_mulmod(t, pk->A[j][i][k], r[j][k], irred[k], ctx_q);
				fmpz_mod_poly_add(c->v[i][k], c->v[i][k], t, ctx_q);
			}
		}
	}

	// Lift m from Rp to Rq. */
	fmpz_mod_poly_get_fmpz_poly(s, m, ctx_p);
	fmpz_set_ui(p2, MODP >> 1);
	for (int i = 0; i < DEGREE; i++) {
		fmpz_poly_get_coeff_fmpz(coeff, s, i);
		if (fmpz_cmp(coeff, p2) >= 0) {
			fmpz_sub(coeff, coeff, p);
		}
		fmpz_mod_poly_set_coeff_fmpz(_m, i, coeff, ctx_q);
	}

	for (int i = 0; i < DIM; i++) {
		fmpz_mod_poly_zero(c->w[i], ctx_q);
		for (int j = 0; j < 2; j++) {
			fmpz_mod_poly_add(c->v[i][j], c->v[i][j], e[i][j], ctx_q);
			fmpz_mod_poly_scalar_mul_fmpz(c->v[i][j], c->v[i][j], p, ctx_q);

			fmpz_mod_poly_mulmod(t, pk->t[j][i], r[j][i], irred[i], ctx_q);
			fmpz_mod_poly_add(c->w[i], c->w[i], t, ctx_q);
		}
		fmpz_mod_poly_add(c->w[i], c->w[i], e_[i], ctx_q);
		fmpz_mod_poly_scalar_mul_fmpz(c->w[i], c->w[i], p, ctx_q);
		fmpz_mod_poly_rem(t, _m, irred[i], ctx_q);
		fmpz_mod_poly_add(c->w[i], c->w[i], t, ctx_q);
	}
	fmpz_mod_poly_clear(_m, ctx_q);
	fmpz_mod_poly_clear(t, ctx_q);
	fmpz_clear(p2);
	fmpz_clear(coeff);
	fmpz_poly_clear(s);
}

// Encrypt a message under a public key.
void encrypt_doit(ciphertext_t *c, fmpz_mod_poly_t m, publickey_t *pk, flint_rand_t rand) {
	qcrt_poly_t r[DIM], e[DIM], e_;

	for (int i = 0; i < DIM; i++) {
		fmpz_mod_poly_init(e_[i], ctx_q);
		for (int j = 0; j < 2; j++) {
			fmpz_mod_poly_init(r[i][j], ctx_q);
			fmpz_mod_poly_init(e[i][j], ctx_q);
		}
		encrypt_sample_short_crt(r[i], ctx_q);
		encrypt_sample_short_crt(e[i], ctx_q);
	}
	encrypt_sample_short_crt(e_, ctx_q);

	encrypt_make(c, r, e, e_, m, pk);

	for (int i = 0; i < DIM; i++) {
		fmpz_mod_poly_clear(e_[i], ctx_q);
		for (int j = 0; j < 2; j++) {
			fmpz_mod_poly_clear(r[i][j], ctx_q);
			fmpz_mod_poly_clear(e[i][j], ctx_q);
		}
	}
}

// Decrypt ciphertext to the original plaintext message.
int encrypt_undo(fmpz_mod_poly_t m, fmpz_mod_poly_t chall, ciphertext_t *c, privatekey_t *sk) {
	fmpz_poly_t s;
	fmpz_mod_poly_t t, u[2];
	int result = 1;

	fmpz_poly_init(s);
	fmpz_mod_poly_init(t, ctx_q);

	for (int i = 0; i < DIM; i++) {
		fmpz_mod_poly_init(u[i], ctx_q);
		fmpz_mod_poly_zero(u[i], ctx_q);
		for (int j = 0; j < 2; j++) {
			fmpz_mod_poly_mulmod(t, c->v[j][i], sk->s1[j][i], irred[i], ctx_q);
			fmpz_mod_poly_add(u[i], u[i], t, ctx_q);
		}
		fmpz_mod_poly_sub(u[i], c->w[i], u[i], ctx_q);
	}
	qcrt_poly_rec(t, u);

	if (chall != NULL) {
		fmpz_mod_poly_mulmod(t, t, chall, large_poly, ctx_q);
	}

	fmpz_mod_poly_get_fmpz_poly(s, t, ctx_q);
	fmpz_t coeff, q2;
	fmpz_init(coeff);
	fmpz_init(q2);

	fmpz_set_str(q2, Q2, 10);
	for (int i = 0; i < DEGREE; i++) {
		fmpz_poly_get_coeff_fmpz(coeff, s, i);
		if (fmpz_cmp(coeff, q2) >= 0) {
			fmpz_sub(coeff, coeff, q);
		}
		fmpz_poly_set_coeff_fmpz(s, i, coeff);
	}

	fmpz_mod_poly_set_fmpz_poly(m, s, ctx_p);

	if (chall != NULL) {
		// Check linf-norm.
		fmpz_set_ui(q2, 12 * SIGMA_E);
		for (int i = 0; i < DEGREE; i++) {
			fmpz_mod_poly_get_coeff_fmpz(coeff, m, i, ctx_p);
			if (fmpz_cmp(coeff, q2) >= 0) {
				result = 0;
			}
		}		
	}

	fmpz_clear(coeff);
	fmpz_poly_clear(s);
	fmpz_mod_poly_clear(t, ctx_q);
	fmpz_mod_poly_clear(u[0], ctx_q);
	fmpz_mod_poly_clear(u[1], ctx_q);
	return result;
}

// Free ciphertext
void encrypt_free(ciphertext_t *c) {
	for (int i = 0; i < DIM; i++) {
		fmpz_mod_poly_clear(c->w[i], ctx_q);
		for (int j = 0; j < DIM; j++) {
			fmpz_mod_poly_clear(c->v[i][j], ctx_q);
		}
	}
}

#ifdef MAIN
// Tests and benchmarks below.
static void test(flint_rand_t rand) {
	publickey_t pk;
	privatekey_t sk;
	ciphertext_t c;
	fmpz_mod_poly_t m, _m, w[2];

	fmpz_mod_poly_init(m, ctx_q);
	fmpz_mod_poly_init(_m, ctx_q);
	fmpz_mod_poly_init(w[0], ctx_q);
	fmpz_mod_poly_init(w[1], ctx_q);

	TEST_BEGIN("CRT representation is correct") {
		fmpz_mod_poly_randtest(m, rand, DEGREE, ctx_q);
		for (int i = 0; i < 2; i++) {
			fmpz_mod_poly_rem(w[i], m, irred[i], ctx_q);
		}
		qcrt_poly_rec(_m, w);
		TEST_ASSERT(fmpz_mod_poly_equal(m, _m, ctx_q) == 1, end);
	} TEST_END;

	fmpz_mod_poly_clear(m, ctx_q);
	fmpz_mod_poly_clear(_m, ctx_q);

	fmpz_mod_poly_init(m, ctx_p);
	fmpz_mod_poly_init(_m, ctx_p);

	TEST_BEGIN("encryption and decryption are consistent") {
		encrypt_sample_short(m, ctx_p);
		encrypt_keygen(&pk, &sk, rand);
		encrypt_doit(&c, m, &pk, rand);
		TEST_ASSERT(encrypt_undo(_m, NULL, &c, &sk) == 1, end);
		TEST_ASSERT(fmpz_mod_poly_equal(m, _m, ctx_p) == 1, end);
	} TEST_END;
end:
	fmpz_mod_poly_clear(w[0], ctx_q);
	fmpz_mod_poly_clear(w[1], ctx_q);
	fmpz_mod_poly_clear(m, ctx_p);
	fmpz_mod_poly_clear(_m, ctx_p);
	encrypt_keyfree(&pk, &sk);
}

static void bench(flint_rand_t rand) {
	publickey_t pk;
	privatekey_t sk;
	ciphertext_t c;
	fmpz_mod_poly_t m, _m;

	fmpz_mod_poly_init(m, ctx_p);
	fmpz_mod_poly_init(_m, ctx_p);

	encrypt_sample_short(m, ctx_p);
	encrypt_keygen(&pk, &sk, rand);

	BENCH_BEGIN("encrypt_doit") {
		BENCH_ADD(encrypt_doit(&c, m, &pk, rand));
	} BENCH_END;

	BENCH_BEGIN("encrypt_undo") {
		BENCH_ADD(encrypt_undo(_m, NULL, &c, &sk));
	} BENCH_END;

	fmpz_mod_poly_clear(m, ctx_p);
	fmpz_mod_poly_clear(_m, ctx_p);
	encrypt_keyfree(&pk, &sk);
}

int main(int argc, char *argv[]) {
	flint_rand_t rand;

	encrypt_setup();

	flint_randinit(rand);

	printf("\n** Tests for lattice-based encryption:\n\n");
	test(rand);

	printf("\n** Benchmarks for lattice-based encryption:\n\n");
	bench(rand);

	encrypt_finish();
}
#endif
