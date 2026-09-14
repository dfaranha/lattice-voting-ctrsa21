/**
 * @file
 *
 * Implementation of the LNP-style proof machinery for the is_bin sub-proof.
 *
 * @ingroup lnp
 */

#include <assert.h>

#include "param.h"
#include "lnp.h"
#include "test.h"
#include "bench.h"
#include "fastrandombytes.h"
#include "sha.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/* The automorphism sigma_k maps X^i to X^(k*i). Exponents live modulo 2*DEGREE
 * because X^DEGREE = -1, so an exponent that lands in [DEGREE, 2*DEGREE) comes
 * back with a sign flip. */
static void auto_exponent(slong *e, int *neg, slong i, slong k) {
	slong t = (i % (2 * DEGREE)) * (k % (2 * DEGREE)) % (2 * DEGREE);

	if (t < 0) {
		t += 2 * DEGREE;
	}
	*neg = (t >= DEGREE);
	*e = *neg ? t - DEGREE : t;
}

/* Absorb a polynomial into a hash by value. FLINT keeps an allocated capacity
 * that is unrelated to the value, so this walks DEGREE coefficients rather
 * than hashing the underlying buffer. */
static void hash_poly(SHA256Context *sha, nmod_poly_t p) {
	uint64_t buf[DEGREE];

	for (int i = 0; i < DEGREE; i++) {
		buf[i] = nmod_poly_get_coeff_ui(p, i);
	}
	SHA256Input(sha, (const uint8_t *)buf, sizeof(buf));
}

/* Derive the Fiat-Shamir challenge from the key, the statement and the whole
 * first message. Slot 3 of the commitment carries the garbage term and belongs
 * to the first message, so it must be absorbed too: leaving it out would let a
 * prover pick it after seeing the challenge. */
static void quad_hash(pcrt_poly_t d, lnpkey_t *key, lnpcom_t *com,
		pcrt_poly_t w[HEIGHT], pcrt_poly_t t) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	uint32_t buf;
	nmod_poly_t c;

	SHA256Reset(&sha);
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < LNP_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				hash_poly(&sha, key->B1[i][j][k]);
			}
		}
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int j = 0; j < LNP_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				hash_poly(&sha, key->b2[i][j][k]);
			}
		}
	}
	for (int k = 0; k < NCRT; k++) {
		for (int i = 0; i < HEIGHT; i++) {
			hash_poly(&sha, com->c1[i][k]);
			hash_poly(&sha, w[i][k]);
		}
		for (int i = 0; i < SLOTS; i++) {
			hash_poly(&sha, com->c2[i][k]);
		}
		hash_poly(&sha, t[k]);
	}
	SHA256Result(&sha, hash);

	/* Sample a challenge with NONZERO coefficients set, as commit.c does, then
	 * reduce it into CRT representation. */
	nmod_poly_init(c, MODP);
	fastrandombytes_setseed(hash);
	nmod_poly_fit_length(c, DEGREE);
	for (int i = 0; i < NONZERO; i++) {
		fastrandombytes((unsigned char *)&buf, sizeof(buf));
		buf = buf % DEGREE;
		while (nmod_poly_get_coeff_ui(c, buf) != 0) {
			fastrandombytes((unsigned char *)&buf, sizeof(buf));
			buf = buf % DEGREE;
		}
		nmod_poly_set_coeff_ui(c, buf, 1);
	}
	pcrt_poly_reduce(d[0], c, 0);
	pcrt_poly_reduce(d[1], c, 1);
	nmod_poly_clear(c);
}

/* Inner product <b, x> over CRT-represented vectors. */
static void inner(pcrt_poly_t out, pcrt_poly_t b[], pcrt_poly_t x[], int len) {
	nmod_poly_t tmp;

	nmod_poly_init(tmp, MODP);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_zero(out[k]);
	}
	for (int i = 0; i < len; i++) {
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(tmp, b[i][k], x[i][k], k);
			nmod_poly_add(out[k], out[k], tmp);
		}
	}
	nmod_poly_clear(tmp);
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void lnp_keyinit(lnpkey_t *key) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < LNP_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(key->B1[i][j][k], MODP);
			}
		}
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int j = 0; j < LNP_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(key->b2[i][j][k], MODP);
			}
		}
	}
}

void lnp_keyfree(lnpkey_t *key) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < LNP_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(key->B1[i][j][k]);
			}
		}
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int j = 0; j < LNP_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(key->b2[i][j][k]);
			}
		}
	}
}

void lnp_keygen(lnpkey_t *key, flint_rand_t rand) {
	/* Hermite normal form: the Ajtai rows start with the identity, and each
	 * message row puts its 1 in the randomness component reserved for it. The
	 * remaining LNP_RANK components carry the MLWE secret that hides all of
	 * them. */
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < LNP_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_zero(key->B1[i][j][k]);
				if (i == j) {
					nmod_poly_set_coeff_ui(key->B1[i][j][k], 0, 1);
				} else if (j >= HEIGHT) {
					commit_sample_rand(key->B1[i][j][k], rand, DEGCRT);
				}
			}
		}
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int j = 0; j < LNP_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_zero(key->b2[i][j][k]);
				if (j == HEIGHT + i) {
					nmod_poly_set_coeff_ui(key->b2[i][j][k], 0, 1);
				} else if (j >= HEIGHT + SLOTS) {
					commit_sample_rand(key->b2[i][j][k], rand, DEGCRT);
				}
			}
		}
	}
}

void lnp_com_init(lnpcom_t *com) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(com->c1[i][k], MODP);
		}
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(com->c2[i][k], MODP);
		}
	}
}

void lnp_com_free(lnpcom_t *com) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(com->c1[i][k]);
		}
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(com->c2[i][k]);
		}
	}
}

void lnp_proof_init(lnpproof_t *pi) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(pi->w[i][k], MODP);
		}
	}
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(pi->z[i][k], MODP);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(pi->t[k], MODP);
	}
}

void lnp_proof_free(lnpproof_t *pi) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(pi->w[i][k]);
		}
	}
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(pi->z[i][k]);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(pi->t[k]);
	}
}

void lnp_commit(lnpcom_t *com, pcrt_poly_t m[SLOTS], lnpkey_t *key,
		pcrt_poly_t r[LNP_WIDTH]) {
	for (int i = 0; i < HEIGHT; i++) {
		inner(com->c1[i], key->B1[i], r, LNP_WIDTH);
	}
	for (int i = 0; i < SLOTS; i++) {
		inner(com->c2[i], key->b2[i], r, LNP_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_add(com->c2[i][k], com->c2[i][k], m[i][k]);
		}
	}
}

int lnp_auto_swaps(slong k) {
	slong r = k % 4;

	if (r < 0) {
		r += 4;
	}
	/* X^(DEGREE/2) squares to -1, so raising it to an odd power k gives
	 * X^(DEGREE/2) when k = 1 mod 4 and -X^(DEGREE/2) when k = 3 mod 4. Since
	 * the two factors are X^(DEGREE/2) + P0 and X^(DEGREE/2) - P0, the second
	 * case sends one factor onto the other. */
	assert(r == 1 || r == 3);
	return (r == 3);
}

void lnp_auto(nmod_poly_t c, const nmod_poly_t a, slong k) {
	nmod_poly_t t;
	slong e;
	int neg;

	assert(k % 2 != 0);

	nmod_poly_init(t, MODP);
	nmod_poly_fit_length(t, DEGREE);
	for (slong i = 0; i < DEGREE; i++) {
		ulong v = nmod_poly_get_coeff_ui(a, i);

		if (v == 0) {
			continue;
		}
		auto_exponent(&e, &neg, i, k);
		if (neg) {
			v = nmod_neg(v, a->mod);
		}
		/* sigma_k is a bijection on exponents, so no two terms collide and a
		 * plain store would do; adding keeps this correct if a ever carries
		 * more than DEGREE coefficients. */
		nmod_poly_set_coeff_ui(t, e,
				nmod_add(nmod_poly_get_coeff_ui(t, e), v, a->mod));
	}
	nmod_poly_set(c, t);
	nmod_poly_clear(t);
}

void lnp_auto_crt(pcrt_poly_t c, pcrt_poly_t a, slong k) {
	nmod_poly_t t[NCRT];
	int swap = lnp_auto_swaps(k);

	/* The residue of a modulo one factor determines the residue of sigma_k(a)
	 * modulo the image factor, because sigma_k maps the one ideal onto the
	 * other. So it suffices to apply sigma_k to the canonical lift of each
	 * component and reduce the image into the target component. */
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(t[i], MODP);
		lnp_auto(t[i], a[i], k);
	}
	for (int i = 0; i < NCRT; i++) {
		int j = swap ? (NCRT - 1 - i) : i;

		pcrt_poly_reduce(c[j], t[i], j);
	}
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(t[i]);
	}
}

/*
 * The quadratic proof. Writing u_j = <b2[j], z> - c * c2[j], an honest
 * transcript has u_j = v_j - c * m_j where v_j = <b2[j], y> is fixed before
 * the challenge. Then
 *
 *   u_0 * u_1 + c * u_2 + u_3
 *     = v_0 * v_1 + v_3 + c * (v_2 - v_0 m_1 - v_1 m_0 - G)
 *                       + c^2 * (m_0 m_1 - m_2),
 *
 * so committing the garbage term G = v_2 - v_0 m_1 - v_1 m_0 in slot 3 kills
 * the linear term, and the quadratic term vanishes exactly when the relation
 * holds. What is left, v_0 v_1 + v_3, is fixed before the challenge and is
 * sent as t, masked by the fresh v_3.
 */
void lnp_quad_prover(lnpproof_t *pi, lnpcom_t *com, pcrt_poly_t m[3],
		lnpkey_t *key, pcrt_poly_t r[LNP_WIDTH]) {
	pcrt_poly_t y[LNP_WIDTH], cr[LNP_WIDTH], v[SLOTS], d, g;
	nmod_poly_t tmp;
	int rej;
	/* sigma_C = 11 * nu * beta * sqrt(k * N), as in the shuffle proof. */
	uint64_t sigma_sqr = 11 * NONZERO * BETA;

	sigma_sqr *= sigma_sqr * DEGREE * LNP_WIDTH;

	nmod_poly_init(tmp, MODP);
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(y[i][k], MODP);
			nmod_poly_init(cr[i][k], MODP);
		}
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(v[i][k], MODP);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(d[k], MODP);
		nmod_poly_init(g[k], MODP);
	}

	do {
		for (int i = 0; i < LNP_WIDTH; i++) {
			commit_sample_gauss_crt(y[i]);
		}
		for (int i = 0; i < HEIGHT; i++) {
			inner(pi->w[i], key->B1[i], y, LNP_WIDTH);
		}
		for (int i = 0; i < SLOTS; i++) {
			inner(v[i], key->b2[i], y, LNP_WIDTH);
		}

		/* G = v_2 - v_0 * m_1 - v_1 * m_0. */
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_set(g[k], v[2][k]);
			pcrt_poly_mulmod(tmp, v[0][k], m[1][k], k);
			nmod_poly_sub(g[k], g[k], tmp);
			pcrt_poly_mulmod(tmp, v[1][k], m[0][k], k);
			nmod_poly_sub(g[k], g[k], tmp);
		}
		/* Commit the garbage term in slot 3, under the same randomness. */
		inner(com->c2[3], key->b2[3], r, LNP_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_add(com->c2[3][k], com->c2[3][k], g[k]);
		}
		/* t = v_0 * v_1 + v_3. */
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(pi->t[k], v[0][k], v[1][k], k);
			nmod_poly_add(pi->t[k], pi->t[k], v[3][k]);
		}

		quad_hash(d, key, com, pi->w, pi->t);

		for (int i = 0; i < LNP_WIDTH; i++) {
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_mulmod(cr[i][k], d[k], r[i][k], k);
				nmod_poly_add(pi->z[i][k], y[i][k], cr[i][k]);
			}
		}
		rej = commit_rej_sampling(pi->z, cr, sigma_sqr, LNP_WIDTH);
	} while (rej);

	nmod_poly_clear(tmp);
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(y[i][k]);
			nmod_poly_clear(cr[i][k]);
		}
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(v[i][k]);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(d[k]);
		nmod_poly_clear(g[k]);
	}
}

int lnp_quad_verifier(lnpproof_t *pi, lnpcom_t *com, lnpkey_t *key) {
	pcrt_poly_t d, u[SLOTS], lhs, rhs;
	nmod_poly_t tmp, rec;
	int result = 1;

	nmod_poly_init(tmp, MODP);
	nmod_poly_init(rec, MODP);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(d[k], MODP);
		nmod_poly_init(lhs[k], MODP);
		nmod_poly_init(rhs[k], MODP);
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(u[i][k], MODP);
		}
	}

	quad_hash(d, key, com, pi->w, pi->t);

	/* The masked opening must be short, or the binding reduction says nothing.
	 * This is a soundness check, so it is an ordinary check and not an
	 * assertion that vanishes under NDEBUG. */
	for (int i = 0; i < LNP_WIDTH; i++) {
		pcrt_poly_rec(rec, pi->z[i]);
		result &= commit_norm2_leq(rec,
				(uint64_t) 4 * DEGREE * SIGMA_C * SIGMA_C);
	}

	/* B1 * z = w + c * c1. */
	for (int i = 0; i < HEIGHT; i++) {
		inner(lhs, key->B1[i], pi->z, LNP_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(tmp, d[k], com->c1[i][k], k);
			nmod_poly_add(rhs[k], pi->w[i][k], tmp);
			result &= nmod_poly_equal(lhs[k], rhs[k]);
		}
	}

	/* u_j = <b2[j], z> - c * c2[j]. */
	for (int i = 0; i < SLOTS; i++) {
		inner(u[i], key->b2[i], pi->z, LNP_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(tmp, d[k], com->c2[i][k], k);
			nmod_poly_sub(u[i][k], u[i][k], tmp);
		}
	}

	/* u_0 * u_1 + c * u_2 + u_3 = t. */
	for (int k = 0; k < NCRT; k++) {
		pcrt_poly_mulmod(lhs[k], u[0][k], u[1][k], k);
		pcrt_poly_mulmod(tmp, d[k], u[2][k], k);
		nmod_poly_add(lhs[k], lhs[k], tmp);
		nmod_poly_add(lhs[k], lhs[k], u[3][k]);
		result &= nmod_poly_equal(lhs[k], pi->t[k]);
	}

	nmod_poly_clear(tmp);
	nmod_poly_clear(rec);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(d[k]);
		nmod_poly_clear(lhs[k]);
		nmod_poly_clear(rhs[k]);
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(u[i][k]);
		}
	}
	return result;
}

/*============================================================================*/
/* Tests                                                                      */
/*============================================================================*/

/* sigma_{-1} in coefficient representation, written out directly from the
 * definition as an independent check on the general routine. */
static void auto_inverse_ref(nmod_poly_t c, nmod_poly_t a) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	nmod_poly_fit_length(t, DEGREE);
	nmod_poly_set_coeff_ui(t, 0, nmod_poly_get_coeff_ui(a, 0));
	for (slong i = 1; i < DEGREE; i++) {
		/* X^-i = -X^(DEGREE-i) because X^DEGREE = -1. */
		nmod_poly_set_coeff_ui(t, DEGREE - i,
				nmod_neg(nmod_poly_get_coeff_ui(a, i), a->mod));
	}
	nmod_poly_set(c, t);
	nmod_poly_clear(t);
}

static void test_auto(flint_rand_t rng) {
	nmod_poly_t a, b, d, e;
	pcrt_poly_t ca, cb;
	/* 2*DEGREE-1 = -1 mod 2*DEGREE, so this is sigma_{-1}. */
	const slong minus1 = 2 * DEGREE - 1;

	nmod_poly_init(a, MODP);
	nmod_poly_init(b, MODP);
	nmod_poly_init(d, MODP);
	nmod_poly_init(e, MODP);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(ca[i], MODP);
		nmod_poly_init(cb[i], MODP);
	}

	TEST_BEGIN("sigma_{-1} matches its definition") {
		commit_sample_rand(a, rng, DEGREE);
		lnp_auto(b, a, minus1);
		auto_inverse_ref(d, a);
		TEST_ASSERT(nmod_poly_equal(b, d) == 1, end);
	} TEST_END;

	TEST_BEGIN("sigma_k is a ring homomorphism") {
		commit_sample_rand(a, rng, DEGREE);
		commit_sample_rand(b, rng, DEGREE);
		/* sigma_k(a*b) = sigma_k(a)*sigma_k(b). */
		commit_poly_mulmod(d, a, b);
		lnp_auto(d, d, minus1);
		lnp_auto(e, a, minus1);
		lnp_auto(b, b, minus1);
		commit_poly_mulmod(e, e, b);
		TEST_ASSERT(nmod_poly_equal(d, e) == 1, end);
	} TEST_END;

	TEST_BEGIN("sigma_{-1} is an involution") {
		commit_sample_rand(a, rng, DEGREE);
		lnp_auto(b, a, minus1);
		lnp_auto(b, b, minus1);
		TEST_ASSERT(nmod_poly_equal(a, b) == 1, end);
	} TEST_END;

	TEST_BEGIN("constant coefficient of sigma_{-1}(a)*b is <a, b>") {
		commit_sample_rand(a, rng, DEGREE);
		commit_sample_rand(b, rng, DEGREE);
		lnp_auto(d, a, minus1);
		commit_poly_mulmod(d, d, b);
		nmod_poly_zero(e);
		for (slong i = 0; i < DEGREE; i++) {
			nmod_poly_set_coeff_ui(e, 0,
					nmod_add(nmod_poly_get_coeff_ui(e, 0),
					nmod_mul(nmod_poly_get_coeff_ui(a, i),
					nmod_poly_get_coeff_ui(b, i), a->mod), a->mod));
		}
		TEST_ASSERT(nmod_poly_get_coeff_ui(d, 0) ==
				nmod_poly_get_coeff_ui(e, 0), end);
	} TEST_END;

	TEST_BEGIN("CRT-aware sigma_k agrees with the coefficient reference") {
		for (slong k = 1; k < 8; k += 2) {
			commit_sample_rand(a, rng, DEGREE);
			/* Reference: map in coefficient representation, then reduce. */
			lnp_auto(b, a, k);
			pcrt_poly_reduce(cb[0], b, 0);
			pcrt_poly_reduce(cb[1], b, 1);
			/* Under test: reduce first, then map in CRT representation. */
			pcrt_poly_reduce(ca[0], a, 0);
			pcrt_poly_reduce(ca[1], a, 1);
			lnp_auto_crt(ca, ca, k);
			TEST_ASSERT(nmod_poly_equal(ca[0], cb[0]) == 1, end);
			TEST_ASSERT(nmod_poly_equal(ca[1], cb[1]) == 1, end);
		}
	} TEST_END;

	TEST_ONCE("sigma_{-1} exchanges the two CRT factors") {
		/* Applying sigma_{-1} to the first irreducible factor must give a
		 * polynomial that vanishes modulo the second, not modulo the first. */
		lnp_auto(a, *commit_irred(0), minus1);
		pcrt_poly_reduce(ca[0], a, 0);
		pcrt_poly_reduce(ca[1], a, 1);
		TEST_ASSERT(nmod_poly_is_zero(ca[1]) == 1, end);
		TEST_ASSERT(nmod_poly_is_zero(ca[0]) == 0, end);
		TEST_ASSERT(lnp_auto_swaps(minus1) == 1, end);
		TEST_ASSERT(lnp_auto_swaps(5) == 0, end);
	} TEST_END;

  end:
	nmod_poly_clear(a);
	nmod_poly_clear(b);
	nmod_poly_clear(d);
	nmod_poly_clear(e);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(ca[i]);
		nmod_poly_clear(cb[i]);
	}
}

static void test_quad(flint_rand_t rng) {
	lnpkey_t key;
	lnpcom_t com;
	lnpproof_t pi;
	pcrt_poly_t r[LNP_WIDTH], m[SLOTS];
	nmod_poly_t a, b;

	lnp_keyinit(&key);
	lnp_com_init(&com);
	lnp_proof_init(&pi);
	lnp_keygen(&key, rng);
	nmod_poly_init(a, MODP);
	nmod_poly_init(b, MODP);
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(r[i][k], MODP);
		}
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(m[i][k], MODP);
		}
	}

	TEST_BEGIN("quadratic proof accepts an honest product") {
		for (int i = 0; i < LNP_WIDTH; i++) {
			commit_sample_short_crt(r[i]);
		}
		commit_sample_rand_crt(m[0], rng);
		commit_sample_rand_crt(m[1], rng);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(m[2][k], m[0][k], m[1][k], k);
			nmod_poly_zero(m[3][k]);
		}
		lnp_commit(&com, m, &key, r);
		lnp_quad_prover(&pi, &com, m, &key, r);
		TEST_ASSERT(lnp_quad_verifier(&pi, &com, &key) == 1, end);
	} TEST_END;

	TEST_BEGIN("quadratic proof rejects a wrong product") {
		for (int i = 0; i < LNP_WIDTH; i++) {
			commit_sample_short_crt(r[i]);
		}
		commit_sample_rand_crt(m[0], rng);
		commit_sample_rand_crt(m[1], rng);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(m[2][k], m[0][k], m[1][k], k);
			/* Claim a product that is off by one. */
			nmod_poly_set_coeff_ui(m[2][k], 0,
					nmod_add(nmod_poly_get_coeff_ui(m[2][k], 0), 1, m[2][k]->mod));
			nmod_poly_zero(m[3][k]);
		}
		lnp_commit(&com, m, &key, r);
		lnp_quad_prover(&pi, &com, m, &key, r);
		TEST_ASSERT(lnp_quad_verifier(&pi, &com, &key) == 0, end);
	} TEST_END;

	TEST_BEGIN("quadratic proof rejects a tampered opening") {
		for (int i = 0; i < LNP_WIDTH; i++) {
			commit_sample_short_crt(r[i]);
		}
		commit_sample_rand_crt(m[0], rng);
		commit_sample_rand_crt(m[1], rng);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(m[2][k], m[0][k], m[1][k], k);
			nmod_poly_zero(m[3][k]);
		}
		lnp_commit(&com, m, &key, r);
		lnp_quad_prover(&pi, &com, m, &key, r);
		nmod_poly_set_coeff_ui(pi.z[0][0], 0,
				nmod_add(nmod_poly_get_coeff_ui(pi.z[0][0], 0), 1,
				pi.z[0][0]->mod));
		TEST_ASSERT(lnp_quad_verifier(&pi, &com, &key) == 0, end);
	} TEST_END;

  end:
	nmod_poly_clear(a);
	nmod_poly_clear(b);
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(r[i][k]);
		}
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(m[i][k]);
		}
	}
	lnp_keyfree(&key);
	lnp_com_free(&com);
	lnp_proof_free(&pi);
}

/* Select which phases to run: "test", "bench", or neither for both. Same
 * helper as in commit.c and shuffle.c; it stays local to each binary so that
 * the test harness does not have to be shared between them. */
static int phase_selected(int argc, char *argv[], const char *phase) {
	return argc < 2 || strcmp(argv[1], phase) == 0;
}

int main(int argc, char *argv[]) {
	flint_rand_t rand;

	flint_rand_init(rand);
	commit_setup();

	if (phase_selected(argc, argv, "test")) {
		printf("\n** Tests for the LNP proof machinery:\n\n");
		test_auto(rand);
		test_quad(rand);
	}

	commit_finish();
	flint_rand_clear(rand);
	return 0;
}
