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
#define Q		"18446744073709551557"
#define Q0		"2296021864060584341"
#define Q1		"16150722209648967216"

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

/* Scratch space for the multiplication routines. Like the rest of this module,
 * they are not reentrant. */
static fmpz_mod_poly_t mul_tmp;

/**
 * Reduce a polynomial modulo (x^m + r) in a given context.
 *
 * Both moduli used here have this shape: the cyclotomic polynomial is
 * x^DEGREE + 1, and each CRT factor is x^DEGCRT + q_i. Reduction is then a fold
 * of each successive block of m coefficients onto the lowest one, scaled by
 * (-r)^k, which avoids the generic polynomial division that
 * fmpz_mod_poly_mulmod performs.
 *
 * The input may have any degree. The output may alias the input: the copy loop
 * is then a no-op, and the fold loop only reads coefficients at or above m,
 * which it never writes.
 *
 * @param[out] c		- the reduced polynomial.
 * @param[in] a			- the polynomial to reduce.
 * @param[in] m			- the degree of the modulus.
 * @param[in] r			- the constant coefficient of the modulus.
 * @param[in] ctx		- the context for modular arithmetic.
 */
static void fold_mod(fmpz_mod_poly_t c, const fmpz_mod_poly_t a, slong m,
		const fmpz_t r, const fmpz_mod_ctx_t ctx) {
	fmpz_t f, t;
	slong len = a->length;

	fmpz_init_set_ui(f, 1);
	fmpz_init(t);

	fmpz_mod_poly_fit_length(c, m, ctx);
	for (slong j = 0; j < m; j++) {
		if (j < len) {
			fmpz_set(c->coeffs + j, a->coeffs + j);
		} else {
			fmpz_zero(c->coeffs + j);
		}
	}
	/* The block of coefficients at x^(k*m) contributes a factor (-r)^k. */
	for (slong base = m; base < len; base += m) {
		fmpz_mod_mul(f, f, r, ctx);
		fmpz_mod_neg(f, f, ctx);
		for (slong j = 0; j < m && base + j < len; j++) {
			fmpz_mod_mul(t, a->coeffs + base + j, f, ctx);
			fmpz_mod_add(c->coeffs + j, c->coeffs + j, t, ctx);
		}
	}
	_fmpz_mod_poly_set_length(c, m);
	_fmpz_mod_poly_normalise(c);

	fmpz_clear(f);
	fmpz_clear(t);
}

/* Constant coefficient of the i-th CRT factor, reduced into a context. */
static void irred_const(fmpz_t r, int i, const fmpz_mod_ctx_t ctx) {
	fmpz_mod_poly_get_coeff_fmpz(r, irred[i], 0, ctx_q);
	fmpz_mod_set_fmpz(r, r, ctx);
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

/* Multiply two polynomials modulo the i-th CRT factor. */
void qcrt_poly_mulmod(fmpz_mod_poly_t c, const fmpz_mod_poly_t a,
		const fmpz_mod_poly_t b, int i, const fmpz_mod_ctx_t ctx) {
	fmpz_t r;

	fmpz_init(r);
	irred_const(r, i, ctx);
	fmpz_mod_poly_mul(mul_tmp, a, b, ctx);
	fold_mod(c, mul_tmp, DEGCRT, r, ctx);
	fmpz_clear(r);
}

/* Reduce a polynomial into the i-th CRT component. */
void qcrt_poly_reduce(fmpz_mod_poly_t c, const fmpz_mod_poly_t a, int i,
		const fmpz_mod_ctx_t ctx) {
	fmpz_t r;

	fmpz_init(r);
	irred_const(r, i, ctx);
	fold_mod(c, a, DEGCRT, r, ctx);
	fmpz_clear(r);
}

/* Multiply two polynomials in the cyclotomic ring (x^DEGREE + 1). */
void encrypt_poly_mulmod(fmpz_mod_poly_t c, const fmpz_mod_poly_t a,
		const fmpz_mod_poly_t b, const fmpz_mod_ctx_t ctx) {
	fmpz_t one;

	fmpz_init_set_ui(one, 1);
	fmpz_mod_poly_mul(mul_tmp, a, b, ctx);
	fold_mod(c, mul_tmp, DEGREE, one, ctx);
	fmpz_clear(one);
}

/* Recover polynomial from CRT representation. */
void qcrt_poly_rec(fmpz_mod_poly_t c, qcrt_poly_t a) {
	fmpz_mod_poly_t t;
	fmpz_t one;

	fmpz_mod_poly_init(t, ctx_q);
	fmpz_init_set_ui(one, 1);

	fmpz_mod_poly_sub(t, a[0], a[1], ctx_q);
	fmpz_mod_poly_mul(t, t, inv[1], ctx_q);
	fmpz_mod_poly_mul(c, t, irred[1], ctx_q);
	fmpz_mod_poly_add(c, c, a[1], ctx_q);
	/* The product of the two CRT factors is the cyclotomic polynomial, so the
	 * reduction is a fold: no need to recompute that product and divide. */
	fold_mod(c, c, DEGREE, one, ctx_q);

	fmpz_clear(one);
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
			fmpz_mod_poly_set_coeff_fmpz(r, (i + j / 2) % DEGREE, coeff, ctx);
		}
	}

	fmpz_clear(coeff);
}

// Sample short element in CRT representation.
void encrypt_sample_short_crt(fmpz_mod_poly_t r[2], fmpz_mod_ctx_t ctx) {
	fmpz_mod_poly_t t;

	fmpz_mod_poly_init(t, ctx);
	encrypt_sample_short(t, ctx);
	for (int i = 0; i < NCRT; i++) {
		qcrt_poly_reduce(r[i], t, i, ctx);
	}

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
	fmpz_mod_poly_init(mul_tmp, ctx_q);
	for (int i = 0; i < NCRT; i++) {
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
	fmpz_mod_poly_clear(mul_tmp, ctx_q);
	for (int i = 0; i < NCRT; i++) {
		fmpz_mod_poly_clear(irred[i], ctx_q);
		fmpz_mod_poly_clear(inv[i], ctx_q);
	}
	fmpz_mod_ctx_clear(ctx_p);
	fmpz_mod_ctx_clear(ctx_q);
	fmpz_clear(p);
	fmpz_clear(q);
}

// Generate a key pair.
// Initialise a ciphertext.
void encrypt_cipher_init(ciphertext_t *c) {
	for (int i = 0; i < NCRT; i++) {
		fmpz_mod_poly_init(c->w[i], ctx_q);
	}
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < NCRT; j++) {
			fmpz_mod_poly_init(c->v[i][j], ctx_q);
		}
	}
}

// Initialise a key pair.
void encrypt_keyinit(publickey_t *pk, privatekey_t *sk) {
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < NCRT; j++) {
			fmpz_mod_poly_init(sk->s1[i][j], ctx_q);
			fmpz_mod_poly_init(sk->s2[i][j], ctx_q);
			fmpz_mod_poly_init(pk->t[i][j], ctx_q);
		}
		for (int j = 0; j < DIM; j++) {
			for (int k = 0; k < NCRT; k++) {
				fmpz_mod_poly_init(pk->A[i][j][k], ctx_q);
			}
		}
	}
}

// Generate a key pair.
void encrypt_keygen(publickey_t *pk, privatekey_t *sk, flint_rand_t rand) {
	fmpz_mod_poly_t t;

	fmpz_mod_poly_init(t, ctx_q);
	for (int i = 0; i < DIM; i++) {
		encrypt_sample_short_crt(sk->s1[i], ctx_q);
		encrypt_sample_short_crt(sk->s2[i], ctx_q);
	}
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < NCRT; j++) {
			fmpz_mod_poly_zero(pk->t[i][j], ctx_q);
		}
		for (int j = 0; j < DIM; j++) {
			for (int k = 0; k < NCRT; k++) {
				fmpz_mod_poly_randtest(pk->A[i][j][k], rand, DEGCRT, ctx_q);
			}
		}
		for (int k = 0; k < NCRT; k++) {
			fmpz_mod_poly_add(pk->t[i][k], pk->t[i][k], sk->s2[i][k], ctx_q);
		}
	}

	// Compute (A, t = As_1 + s_2).
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			for (int k = 0; k < NCRT; k++) {
				qcrt_poly_mulmod(t, pk->A[i][j][k], sk->s1[j][k], k, ctx_q);
				fmpz_mod_poly_add(pk->t[i][k], pk->t[i][k], t, ctx_q);
			}
		}
	}
	fmpz_mod_poly_clear(t, ctx_q);
}

// Free key pair.
void encrypt_keyfree(publickey_t *pk, privatekey_t *sk) {
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < NCRT; j++) {
			fmpz_mod_poly_clear(pk->t[i][j], ctx_q);
			fmpz_mod_poly_clear(sk->s1[i][j], ctx_q);
			fmpz_mod_poly_clear(sk->s2[i][j], ctx_q);
		}
		for (int j = 0; j < DIM; j++) {
			for (int k = 0; k < NCRT; k++) {
				fmpz_mod_poly_clear(pk->A[i][j][k], ctx_q);
			}
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
		for (int j = 0; j < NCRT; j++) {
			fmpz_mod_poly_zero(c->v[i][j], ctx_q);
		}
	}

	fmpz_mod_poly_init(t, ctx_q);
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			for (int k = 0; k < NCRT; k++) {
				qcrt_poly_mulmod(t, pk->A[j][i][k], r[j][k], k, ctx_q);
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

	/* v is indexed by the module dimension and then the CRT component. */
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < NCRT; j++) {
			fmpz_mod_poly_add(c->v[i][j], c->v[i][j], e[i][j], ctx_q);
			fmpz_mod_poly_scalar_mul_fmpz(c->v[i][j], c->v[i][j], p, ctx_q);
		}
	}

	/* w has one entry per CRT component, each a sum over the module
	 * dimension. */
	for (int k = 0; k < NCRT; k++) {
		fmpz_mod_poly_zero(c->w[k], ctx_q);
		for (int j = 0; j < DIM; j++) {
			qcrt_poly_mulmod(t, pk->t[j][k], r[j][k], k, ctx_q);
			fmpz_mod_poly_add(c->w[k], c->w[k], t, ctx_q);
		}
		fmpz_mod_poly_add(c->w[k], c->w[k], e_[k], ctx_q);
		fmpz_mod_poly_scalar_mul_fmpz(c->w[k], c->w[k], p, ctx_q);
		qcrt_poly_reduce(t, _m, k, ctx_q);
		fmpz_mod_poly_add(c->w[k], c->w[k], t, ctx_q);
	}
	fmpz_mod_poly_clear(_m, ctx_q);
	fmpz_mod_poly_clear(t, ctx_q);
	fmpz_clear(p2);
	fmpz_clear(coeff);
	fmpz_poly_clear(s);
}

// Encrypt a message under a public key.
void encrypt_doit(ciphertext_t *c, fmpz_mod_poly_t m, publickey_t *pk,
		flint_rand_t rand) {
	qcrt_poly_t r[DIM], e[DIM], e_;

	for (int i = 0; i < NCRT; i++) {
		fmpz_mod_poly_init(e_[i], ctx_q);
	}
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < NCRT; j++) {
			fmpz_mod_poly_init(r[i][j], ctx_q);
			fmpz_mod_poly_init(e[i][j], ctx_q);
		}
		encrypt_sample_short_crt(r[i], ctx_q);
		encrypt_sample_short_crt(e[i], ctx_q);
	}
	encrypt_sample_short_crt(e_, ctx_q);

	encrypt_make(c, r, e, e_, m, pk);

	for (int i = 0; i < NCRT; i++) {
		fmpz_mod_poly_clear(e_[i], ctx_q);
	}
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < NCRT; j++) {
			fmpz_mod_poly_clear(r[i][j], ctx_q);
			fmpz_mod_poly_clear(e[i][j], ctx_q);
		}
	}
}

// Decrypt ciphertext to the original plaintext message.
int encrypt_undo(fmpz_mod_poly_t m, fmpz_mod_poly_t chall, ciphertext_t *c,
		privatekey_t *sk) {
	fmpz_poly_t s;
	fmpz_mod_poly_t t, _t, u[2];
	fmpz_t coeff, q2;
	int result = 1;
	fmpz_init(coeff);
	fmpz_init(q2);

	fmpz_poly_init(s);
	fmpz_mod_poly_init(t, ctx_q);
	fmpz_mod_poly_init(_t, ctx_q);

	for (int i = 0; i < NCRT; i++) {
		fmpz_mod_poly_init(u[i], ctx_q);
		fmpz_mod_poly_zero(u[i], ctx_q);
		for (int j = 0; j < DIM; j++) {
			qcrt_poly_mulmod(t, c->v[j][i], sk->s1[j][i], i, ctx_q);
			fmpz_mod_poly_add(u[i], u[i], t, ctx_q);
		}
		fmpz_mod_poly_sub(u[i], c->w[i], u[i], ctx_q);
	}
	qcrt_poly_rec(t, u);

	if (chall != NULL) {
		fmpz_set_ui(q2, MODP / 2);
		for (int i = 0; i < DEGREE; i++) {
			fmpz_mod_poly_get_coeff_fmpz(coeff, chall, i, ctx_p);
			if (fmpz_cmp(coeff, q2) >= 0) {
				fmpz_sub_ui(coeff, coeff, MODP);
			}
			fmpz_mod_poly_set_coeff_fmpz(_t, i, coeff, ctx_q);
		}
		encrypt_poly_mulmod(t, t, _t, ctx_q);
	}

	fmpz_mod_poly_get_fmpz_poly(s, t, ctx_q);

	/* Centre the coefficients around zero: q2 = floor(q / 2). Deriving this
	 * from q rather than hard-coding it keeps it correct when q changes. */
	fmpz_fdiv_q_ui(q2, q, 2);
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
				//TODO: fixme
				//result = 0;
			}
		}
	}

	fmpz_clear(coeff);
	fmpz_poly_clear(s);
	fmpz_mod_poly_clear(t, ctx_q);
	fmpz_mod_poly_clear(_t, ctx_q);
	fmpz_mod_poly_clear(u[0], ctx_q);
	fmpz_mod_poly_clear(u[1], ctx_q);
	return result;
}

// Free ciphertext
void encrypt_free(ciphertext_t *c) {
	for (int i = 0; i < NCRT; i++) {
		fmpz_mod_poly_clear(c->w[i], ctx_q);
	}
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < NCRT; j++) {
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
		for (int i = 0; i < NCRT; i++) {
			qcrt_poly_reduce(w[i], m, i, ctx_q);
		}
		qcrt_poly_rec(_m, w);
		TEST_ASSERT(fmpz_mod_poly_equal(m, _m, ctx_q) == 1, end);
	} TEST_END;

	fmpz_mod_poly_clear(m, ctx_q);
	fmpz_mod_poly_clear(_m, ctx_q);

	fmpz_mod_poly_init(m, ctx_p);
	fmpz_mod_poly_init(_m, ctx_p);

	encrypt_keyinit(&pk, &sk);
	encrypt_cipher_init(&c);

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
	encrypt_free(&c);
}

static void bench(flint_rand_t rand) {
	publickey_t pk;
	privatekey_t sk;
	ciphertext_t c;
	fmpz_mod_poly_t m, _m;

	fmpz_mod_poly_init(m, ctx_p);
	fmpz_mod_poly_init(_m, ctx_p);

	encrypt_sample_short(m, ctx_p);
	encrypt_keyinit(&pk, &sk);
	encrypt_keygen(&pk, &sk, rand);
	encrypt_cipher_init(&c);

	BENCH_BEGIN("encrypt_doit") {
		BENCH_ADD(encrypt_doit(&c, m, &pk, rand));
	} BENCH_END;

	BENCH_BEGIN("encrypt_undo") {
		BENCH_ADD(encrypt_undo(_m, NULL, &c, &sk));
	} BENCH_END;

	encrypt_free(&c);

	fmpz_mod_poly_clear(m, ctx_p);
	fmpz_mod_poly_clear(_m, ctx_p);
	encrypt_keyfree(&pk, &sk);
}

/* Select which phases to run: "test", "bench", or neither for both. Keeping
 * the benchmarks out of a test run matters in practice, since they dominate
 * the runtime by two orders of magnitude. */
static int phase_selected(int argc, char *argv[], const char *phase) {
	return argc < 2 || strcmp(argv[1], phase) == 0;
}

int main(int argc, char *argv[]) {
	flint_rand_t rand;

	encrypt_setup();

	flint_rand_init(rand);

	if (phase_selected(argc, argv, "test")) {
		printf("\n** Tests for lattice-based encryption:\n\n");
		test(rand);
	}

	if (phase_selected(argc, argv, "bench")) {
		printf("\n** Benchmarks for lattice-based encryption:\n\n");
		bench(rand);
	}

	encrypt_finish();
}
#endif
