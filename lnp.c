/**
 * @file
 *
 * Implementation of the LNP-style proof machinery for the is_bin sub-proof.
 *
 * @ingroup lnp
 */

#include <assert.h>
#include <math.h>
#include <string.h>

#include "param.h"
#include "lnp.h"
#include "test.h"
#include "bench.h"
#include "fastrandombytes.h"
#include "sha.h"
#include "gaussian.h"


/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/


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


/*============================================================================*/
/* Approximate range proof                                                    */
/*============================================================================*/

/* Centre a residue around zero. */
static int64_t centre(ulong a) {
	return (a > MODP / 2) ? (int64_t) a - (int64_t) MODP : (int64_t) a;
}

/* Seed the projection matrix from the commitment. Every row is regenerated on
 * demand from this seed rather than stored, since the matrix is PROJ by DEGREE
 * and is needed twice. */
static void proj_seed(uint8_t hash[SHA256HashSize], lnpkey_t *key,
		lnpcom_t *com, lnpmaskcom_t *mcom) {
	SHA256Context sha;

	SHA256Reset(&sha);
	for (int i = 0; i < SLOTS; i++) {
		for (int j = 0; j < LNP_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				hash_poly(&sha, key->b2[i][j][k]);
			}
		}
	}
	/* Only the Ajtai part and the slots the projection speaks about: the
	 * witness, the projection mask and the constant-coefficient masks. The
	 * garbage slots are deliberately left out. They are written by the product
	 * relation, which depends on the masking, so binding the projection to
	 * them would put the whole projection inside the rejection loop for no
	 * soundness benefit: what the seed has to guarantee is that the prover
	 * cannot predict the matrix before committing the witness and the mask,
	 * and those slots are covered. */
	for (int k = 0; k < NCRT; k++) {
		for (int i = 0; i < HEIGHT; i++) {
			hash_poly(&sha, com->c1[i][k]);
		}
		hash_poly(&sha, com->c2[SLOT_S][k]);
		hash_poly(&sha, com->c2[SLOT_W][k]);
		/* The batch masks live in their own commitment now, and the scalars
		 * derived from this seed are exactly what they have to be fixed
		 * before, so it has to bind them. */
		for (int i = 0; i < HEIGHT; i++) {
			hash_poly(&sha, mcom->c1[i][k]);
		}
		for (int i = 0; i < LNP_LAMBDA; i++) {
			hash_poly(&sha, mcom->c2[i][k]);
		}
	}
	SHA256Result(&sha, hash);
}

/* Draw the next row of the projection matrix: DEGREE signs, one bit each. */
static void proj_row(int8_t row[DEGREE]) {
	uint8_t buf[DEGREE / 8];

	fastrandombytes(buf, sizeof(buf));
	for (int i = 0; i < DEGREE; i++) {
		row[i] = ((buf[i >> 3] >> (i & 7)) & 1) ? 1 : -1;
	}
}

/* Derive the aggregation scalars. These come after the projection has been
 * published, so that the masks committed earlier cannot be adapted to them. */
static void proj_scalars(ulong mu[LNP_LAMBDA][PROJ], ulong nu[LNP_LAMBDA],
		uint8_t seed[SHA256HashSize], ulong z[PROJ]) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	uint64_t buf;

	SHA256Reset(&sha);
	SHA256Input(&sha, seed, SHA256HashSize);
	SHA256Input(&sha, (const uint8_t *)z, PROJ * sizeof(ulong));
	SHA256Result(&sha, hash);

	fastrandombytes_setseed(hash);
	for (int j = 0; j < LNP_LAMBDA; j++) {
		for (int i = 0; i < PROJ; i++) {
			fastrandombytes((unsigned char *)&buf, sizeof(buf));
			mu[j][i] = buf % MODP;
		}
		/* One more scalar per mask, weighting the constant-coefficient claim
		 * that used to be a separate proof. Folding it in here is what lets a
		 * single set of masks cover both, which is also what stops the two
		 * from sharing a mask and cancelling it. */
		fastrandombytes((unsigned char *)&buf, sizeof(buf));
		nu[j] = buf % MODP;
	}
}

/*
 * One pass over the projection matrix.
 *
 * If sc is not NULL, rs[i] is set to <r_i, s>. In every case acc[j] gathers
 * the row combination sum_i mu[j][i] * r_i, whose automorphism is the public
 * multiplier that the consistency relation applies to the committed witness:
 * since sigma is linear, sum_i mu[j][i] * sigma(r_i) is sigma of that sum, so
 * one automorphism per mask suffices instead of one per row.
 */
static void proj_pass(uint8_t seed[SHA256HashSize], ulong *raw,
		ulong mu[LNP_LAMBDA][PROJ], nmod_poly_t sc, ulong *rs) {
	int8_t row[DEGREE];
	nmod_t mod;

	nmod_init(&mod, MODP);

	fastrandombytes_setseed(seed);
	if (raw != NULL) {
		memset(raw, 0, (size_t) LNP_LAMBDA * DEGREE * sizeof(ulong));
	}
	for (int i = 0; i < PROJ; i++) {
		proj_row(row);
		if (sc != NULL) {
			/* Centred coefficients are below 2^39 and there are DEGREE of
			 * them, so this signed accumulation cannot overflow. */
			int64_t dot = 0;

			for (int t = 0; t < DEGREE; t++) {
				dot += row[t] * centre(nmod_poly_get_coeff_ui(sc, t));
			}
			rs[i] = (dot < 0) ? MODP - ((-dot) % MODP) : (ulong)(dot % MODP);
			if (rs[i] == MODP) {
				rs[i] = 0;
			}
		}
		if (raw != NULL && mu != NULL) {
			for (int j = 0; j < LNP_LAMBDA; j++) {
				ulong m = mu[j][i];
				ulong neg = nmod_neg(m, mod);
				ulong *dst = raw + (size_t) j * DEGREE;

				/* Accumulating into a plain array rather than through
				 * nmod_poly_set_coeff_ui matters here: this loop runs
				 * PROJ * LNP_LAMBDA * DEGREE times, and the accessor
				 * renormalises the polynomial on every call. */
				for (int t = 0; t < DEGREE; t++) {
					dst[t] = nmod_add(dst[t], row[t] > 0 ? m : neg, mod);
				}
			}
		}
	}
}

/* Build the public multipliers of the consistency relation: P applies to the
 * witness and sigma(M) to the packed projection mask. */
static void proj_public(pcrt_poly_t P[LNP_LAMBDA], pcrt_poly_t M[LNP_LAMBDA],
		ulong Z[LNP_LAMBDA], uint8_t seed[SHA256HashSize],
		ulong mu[LNP_LAMBDA][PROJ], ulong z[PROJ]) {
	nmod_poly_t acc[LNP_LAMBDA], t;
	const slong minus1 = 2 * DEGREE - 1;
	ulong *raw = (ulong *) flint_malloc((size_t) LNP_LAMBDA * DEGREE *
			sizeof(ulong));

	nmod_poly_init(t, MODP);
	for (int j = 0; j < LNP_LAMBDA; j++) {
		nmod_poly_init(acc[j], MODP);
	}
	proj_pass(seed, raw, mu, NULL, NULL);
	for (int j = 0; j < LNP_LAMBDA; j++) {
		nmod_poly_zero(acc[j]);
		nmod_poly_fit_length(acc[j], DEGREE);
		for (int t2 = 0; t2 < DEGREE; t2++) {
			nmod_poly_set_coeff_ui(acc[j], t2, raw[(size_t) j * DEGREE + t2]);
		}
		/* P_j = sigma(sum_i mu_ji r_i). */
		lnp_auto(acc[j], acc[j], minus1);
		pcrt_poly_reduce(P[j][0], acc[j], 0);
		pcrt_poly_reduce(P[j][1], acc[j], 1);
		/* M_j = sigma(sum_i mu_ji X^i), so that its constant coefficient
		 * against the packed mask picks out sum_i mu_ji w_i. */
		nmod_poly_zero(t);
		nmod_poly_fit_length(t, DEGREE);
		for (int i = 0; i < PROJ; i++) {
			nmod_poly_set_coeff_ui(t, i, mu[j][i]);
		}
		lnp_auto(t, t, minus1);
		pcrt_poly_reduce(M[j][0], t, 0);
		pcrt_poly_reduce(M[j][1], t, 1);
		/* Z_j = sum_i mu_ji z_i. */
		Z[j] = 0;
		for (int i = 0; i < PROJ; i++) {
			Z[j] = nmod_add(Z[j], nmod_mul(mu[j][i], z[i], t->mod), t->mod);
		}
	}
	flint_free(raw);
	nmod_poly_clear(t);
	for (int j = 0; j < LNP_LAMBDA; j++) {
		nmod_poly_clear(acc[j]);
	}
}



void lnp_sample_proj_mask(pcrt_poly_t w, nmod_poly_t raw) {
	nmod_poly_zero(raw);
	nmod_poly_fit_length(raw, DEGREE);
	for (int i = 0; i < PROJ; i++) {
		int64_t c = discrete_gaussian_proj(0.0);

		nmod_poly_set_coeff_ui(raw, i,
				(c < 0) ? MODP - ((-c) % MODP) : (ulong)(c % MODP));
	}
	pcrt_poly_reduce(w[0], raw, 0);
	pcrt_poly_reduce(w[1], raw, 1);
}

/* Squared 2-norm of a centred vector against a bound, without overflowing.
 * Comparing a coefficient against the bound is not enough, since the bound is
 * itself a squared norm and a coefficient below it can still square past
 * 2^64; dividing decides the same question safely. */
static int proj_norm2_leq(ulong z[PROJ], uint64_t bound) {
	uint64_t norm = 0;

	for (int i = 0; i < PROJ; i++) {
		int64_t c = centre(z[i]);
		uint64_t a = (uint64_t)(c < 0 ? -c : c);

		if (a != 0 && a > bound / a) {
			return 0;
		}
		norm += a * a;
		if (norm > bound) {
			return 0;
		}
	}
	return 1;
}

/* Rejection sampling on the published projection.
 *
 * A witness far outside the range makes the projection enormous, and then the
 * products below would overflow. Such a transcript can never be accepted, so
 * it is rejected before the arithmetic rather than after. */
static int proj_reject(ulong z[PROJ], ulong v[PROJ]) {
	double r, u, M = 3.8;
	int64_t dot = 0, norm = 0;
	/* Beyond this the projection cannot pass the verifier's check anyway, and
	 * the accumulation below would no longer fit. */
	const int64_t lim = (int64_t) 1 << 28;

	for (int i = 0; i < PROJ; i++) {
		int64_t c0 = centre(z[i]), c1 = centre(v[i]);

		if (c0 > lim || c0 < -lim || c1 > lim || c1 < -lim) {
			return 1;
		}
		dot += c0 * c1;
		norm += c1 * c1;
	}
	u = commit_uniform_double();
	r = exp((-2.0 * dot + norm) / (2.0 * (double) SIGMA_P * SIGMA_P)) / M;
	return u > r;
}


/* The combined key row and commitment of the j-th consistency relation: P_j
 * applied to the witness slot, M_j to the packed mask slot, and nu_j to the
 * product slot. Both P_j and M_j already carry the automorphism, applied once
 * when proj_public built them. */
static void range_row(pcrt_poly_t b[LNP_WIDTH], pcrt_poly_t t, lnpkey_t *key,
		lnpcom_t *com, pcrt_poly_t P, pcrt_poly_t M, ulong nu) {
	nmod_poly_t tmp;

	nmod_poly_init(tmp, MODP);
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(b[i][k], P[k], key->b2[SLOT_S][i][k], k);
			pcrt_poly_mulmod(tmp, M[k], key->b2[SLOT_W][i][k], k);
			nmod_poly_add(b[i][k], b[i][k], tmp);
			nmod_poly_scalar_mul_nmod(tmp, key->b2[SLOT_F][i][k], nu);
			nmod_poly_add(b[i][k], b[i][k], tmp);
		}
	}
	if (t != NULL) {
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(t[k], P[k], com->c2[SLOT_S][k], k);
			pcrt_poly_mulmod(tmp, M[k], com->c2[SLOT_W][k], k);
			nmod_poly_add(t[k], t[k], tmp);
			nmod_poly_scalar_mul_nmod(tmp, com->c2[SLOT_F][k], nu);
			nmod_poly_add(t[k], t[k], tmp);
		}
	}
	nmod_poly_clear(tmp);
}

/* Reduce a scalar into both CRT components. */
static void scalar_crt(pcrt_poly_t out, ulong a) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	nmod_poly_set_coeff_ui(t, 0, a);
	pcrt_poly_reduce(out[0], t, 0);
	pcrt_poly_reduce(out[1], t, 1);
	nmod_poly_clear(t);
}

/*
 * The range proof. The witness is projected onto PROJ coordinates with a
 * public sign matrix and the projection published masked, so the verifier
 * learns a bound on ||R s|| and hence, by the projection lemma, on ||s||.
 *
 * Tying the projection back to the commitment is what the constant-coefficient
 * machinery is for. Since <r_i, s> is the constant coefficient of
 * sigma(r_i) * s, and the i-th coefficient of the packed mask W is the
 * constant coefficient of sigma(X^i) * W, the whole system of PROJ equations
 * aggregates into LNP_LAMBDA statements
 *
 *    ct(P_j * s + sigma(M_j) * W - Z_j) = 0
 *
 * with P_j, M_j and Z_j public. That is the shape the constant-coefficient
 * proves, with the scalar multipliers replaced by ring ones.
 */
/* Exposed so that a test can play the part of a prover that reads the scalars
 * and then tries to adapt its masks to them. */
void lnp_scalars_for_test(ulong nu[LNP_LAMBDA], lnpkey_t *key, lnpcom_t *com,
		lnpmaskcom_t *mcom, ulong z[PROJ]) {
	ulong (*mu)[PROJ] = (ulong (*)[PROJ]) flint_malloc((size_t) LNP_LAMBDA *
			PROJ * sizeof(ulong));
	uint8_t seed[SHA256HashSize];

	proj_seed(seed, key, com, mcom);
	proj_scalars(mu, nu, seed, z);
	flint_free(mu);
}



/* The value an honest prover puts in slot SLOT_F: sigma_{-1}(s) * (s - ones),
 * whose constant coefficient is the sum over j of s_j (s_j - 1). */
void lnp_isbin_product(pcrt_poly_t f, pcrt_poly_t s) {
	pcrt_poly_t one, ss;
	nmod_poly_t tmp;

	nmod_poly_init(tmp, MODP);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(one[k], MODP);
		nmod_poly_init(ss[k], MODP);
	}
	lnp_ones(one);
	lnp_auto_crt(ss, s, 2 * DEGREE - 1);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_sub(tmp, s[k], one[k]);
		pcrt_poly_mulmod(f[k], ss[k], tmp, k);
	}
	nmod_poly_clear(tmp);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(one[k]);
		nmod_poly_clear(ss[k]);
	}
}

void lnp_ones(pcrt_poly_t out) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	nmod_poly_fit_length(t, DEGREE);
	for (slong i = 0; i < DEGREE; i++) {
		nmod_poly_set_coeff_ui(t, i, 1);
	}
	pcrt_poly_reduce(out[0], t, 0);
	pcrt_poly_reduce(out[1], t, 1);
	nmod_poly_clear(t);
}






void lnp_sample_ct_zero(pcrt_poly_t g, flint_rand_t rand) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	commit_sample_rand(t, rand, DEGREE);
	nmod_poly_set_coeff_ui(t, 0, 0);
	pcrt_poly_reduce(g[0], t, 0);
	pcrt_poly_reduce(g[1], t, 1);
	nmod_poly_clear(t);
}



/* Derive the LNP_LAMBDA aggregation scalars from the commitment alone. They
 * must not depend on anything chosen after the masks are committed, which is
 * what rules out a prover that picks a mask cancelling a non-zero constant
 * coefficient. */








void lnp_binctx_init(lnpbinctx_t *ctx) {
	for (int j = 0; j < LNP_LAMBDA; j++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(ctx->P[j][k], MODP);
			nmod_poly_init(ctx->M[j][k], MODP);
		}
	}
}

void lnp_binctx_free(lnpbinctx_t *ctx) {
	for (int j = 0; j < LNP_LAMBDA; j++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(ctx->P[j][k]);
			nmod_poly_clear(ctx->M[j][k]);
		}
	}
}

void lnp_binproof_init(lnpbinproof_t *pi) {
	(void) pi;
}

void lnp_binproof_free(lnpbinproof_t *pi) {
	(void) pi;
}

void lnp_bin_public(lnpbinctx_t *ctx, lnpbinproof_t *pi, lnpcom_t *com,
		lnpmaskcom_t *mcom, lnpkey_t *key) {
	uint8_t seed[SHA256HashSize];
	ulong (*mu)[PROJ] = (ulong (*)[PROJ]) flint_malloc((size_t) LNP_LAMBDA *
			PROJ * sizeof(ulong));

	proj_seed(seed, key, com, mcom);
	proj_scalars(mu, ctx->nu, seed, pi->zp);
	proj_public(ctx->P, ctx->M, ctx->Z, seed, mu, pi->zp);
	flint_free(mu);
}

int lnp_bin_setup(lnpbinproof_t *pi, lnpbinctx_t *ctx, lnpbatch_t *batch,
		lnpcom_t *com, lnpmaskcom_t *mcom, pcrt_poly_t s, pcrt_poly_t f,
		nmod_poly_t w, lnpkey_t *key) {
	nmod_poly_t sc, tmp;
	pcrt_poly_t zc, wc;
	uint8_t seed[SHA256HashSize];
	ulong *rs;
	int ok = 1;

	rs = (ulong *) flint_malloc(PROJ * sizeof(ulong));
	nmod_poly_init(sc, MODP);
	nmod_poly_init(tmp, MODP);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(zc[k], MODP);
		nmod_poly_init(wc[k], MODP);
	}
	pcrt_poly_rec(sc, s);
	pcrt_poly_reduce(wc[0], w, 0);
	pcrt_poly_reduce(wc[1], w, 1);

	proj_seed(seed, key, com, mcom);
	proj_pass(seed, NULL, NULL, sc, rs);
	for (int i = 0; i < PROJ; i++) {
		pi->zp[i] = nmod_add(nmod_poly_get_coeff_ui(w, i), rs[i], sc->mod);
	}
	if (proj_reject(pi->zp, rs)) {
		/* Return before touching the batch: the caller will recommit and call
		 * again, and a rejected attempt must not leave its share in the
		 * aggregated values. */
		ok = 0;
		goto done;
	}
	lnp_bin_public(ctx, pi, com, mcom, key);

	/* This message's share of h_j = g_j + sum_l (nu_lj f_l + P_lj s_l +
	 * sigma(M_lj) W_l - Z_lj). The mask g_j is added once, by the caller,
	 * since it is shared across the batch. */
	for (int j = 0; j < LNP_LAMBDA; j++) {
		scalar_crt(zc, ctx->Z[j]);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(tmp, ctx->P[j][k], s[k], k);
			nmod_poly_add(batch->h[j][k], batch->h[j][k], tmp);
			pcrt_poly_mulmod(tmp, ctx->M[j][k], wc[k], k);
			nmod_poly_add(batch->h[j][k], batch->h[j][k], tmp);
			nmod_poly_scalar_mul_nmod(tmp, f[k], ctx->nu[j]);
			nmod_poly_add(batch->h[j][k], batch->h[j][k], tmp);
			nmod_poly_sub(batch->h[j][k], batch->h[j][k], zc[k]);
		}
	}

  done:
	flint_free(rs);
	nmod_poly_clear(sc);
	nmod_poly_clear(tmp);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(zc[k]);
		nmod_poly_clear(wc[k]);
	}
	return ok;
}

/* This message's mask-dependent half. The garbage terms are no longer
 * committed per message: each message contributes its rho-weighted share to
 * the two batch-wide accumulators instead, and its share of the batched
 * challenge-free term T. The second garbage term is weighted by sigma(rho)
 * rather than rho, because the verifier's identity applies sigma to that row,
 * and sigma(sigma(rho) g2) is rho sigma(g2). */
void lnp_bin_first(lnpbinctx_t *ctx, lnpbatch_t *batch, lnpcom_t *com,
		pcrt_poly_t s, lnpkey_t *key, pcrt_poly_t r[LNP_WIDTH],
		pcrt_poly_t y[LNP_WIDTH], pcrt_poly_t rho, pcrt_poly_t garb[2]) {
	pcrt_poly_t v0, v1, b[LNP_WIDTH], g1, g2, sv0, ss, one, part, srho;
	nmod_poly_t tmp;
	const slong minus1 = 2 * DEGREE - 1;

	(void) r;
	nmod_poly_init(tmp, MODP);
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(b[i][k], MODP);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(v0[k], MODP);
		nmod_poly_init(v1[k], MODP);
		nmod_poly_init(g1[k], MODP);
		nmod_poly_init(g2[k], MODP);
		nmod_poly_init(sv0[k], MODP);
		nmod_poly_init(ss[k], MODP);
		nmod_poly_init(one[k], MODP);
		nmod_poly_init(part[k], MODP);
		nmod_poly_init(srho[k], MODP);
	}
	lnp_ones(one);
	lnp_auto_crt(ss, s, minus1);
	lnp_auto_crt(srho, rho, minus1);

	inner(v0, key->b2[SLOT_S], y, LNP_WIDTH);
	inner(v1, key->b2[SLOT_F], y, LNP_WIDTH);
	lnp_auto_crt(sv0, v0, minus1);

	/* g1 = -sigma(v0) (s - 1), g2 = sigma(v1 - v0 sigma(s)). */
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_sub(tmp, s[k], one[k]);
		pcrt_poly_mulmod(g1[k], sv0[k], tmp, k);
		nmod_poly_neg(g1[k], g1[k]);
		pcrt_poly_mulmod(tmp, v0[k], ss[k], k);
		nmod_poly_sub(g2[k], v1[k], tmp);
	}
	lnp_auto_crt(g2, g2, minus1);
	for (int k = 0; k < NCRT; k++) {
		pcrt_poly_mulmod(tmp, rho[k], g1[k], k);
		nmod_poly_add(garb[0][k], garb[0][k], tmp);
		pcrt_poly_mulmod(tmp, srho[k], g2[k], k);
		nmod_poly_add(garb[1][k], garb[1][k], tmp);
		/* T gets rho * sigma(v0) * v0. */
		pcrt_poly_mulmod(tmp, sv0[k], v0[k], k);
		pcrt_poly_mulmod(tmp, rho[k], tmp, k);
		nmod_poly_add(batch->T[k], batch->T[k], tmp);
	}

	for (int j = 0; j < LNP_LAMBDA; j++) {
		range_row(b, NULL, key, com, ctx->P[j], ctx->M[j], ctx->nu[j]);
		inner(part, b, y, LNP_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_add(batch->v[j][k], batch->v[j][k], part[k]);
		}
	}

	nmod_poly_clear(tmp);
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(b[i][k]);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(v0[k]);
		nmod_poly_clear(v1[k]);
		nmod_poly_clear(g1[k]);
		nmod_poly_clear(g2[k]);
		nmod_poly_clear(sv0[k]);
		nmod_poly_clear(ss[k]);
		nmod_poly_clear(one[k]);
		nmod_poly_clear(part[k]);
		nmod_poly_clear(srho[k]);
	}
}

/* This message's share of the two aggregated relations. The quadratic part is
 * no longer checked here: with the garbage terms batched, only the sum over
 * messages closes, so each message adds its rho-weighted quadratic terms to
 * batch->acc and lnp_batch_check compares the total against T. */
int lnp_bin_check(lnpbinproof_t *pi, lnpbinctx_t *ctx, lnpbatch_t *batch,
		lnpcom_t *com, lnpkey_t *key, pcrt_poly_t d, pcrt_poly_t z[LNP_WIDTH],
		pcrt_poly_t acc[LNP_LAMBDA], pcrt_poly_t rho) {
	pcrt_poly_t sd, u[SLOTS], su, lhs, one, t, zc, b[LNP_WIDTH];
	nmod_poly_t tmp;
	int result = 1;
	const slong minus1 = 2 * DEGREE - 1;

	nmod_poly_init(tmp, MODP);
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(b[i][k], MODP);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(sd[k], MODP);
		nmod_poly_init(su[k], MODP);
		nmod_poly_init(lhs[k], MODP);
		nmod_poly_init(one[k], MODP);
		nmod_poly_init(t[k], MODP);
		nmod_poly_init(zc[k], MODP);
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(u[i][k], MODP);
		}
	}
	lnp_ones(one);
	lnp_auto_crt(sd, d, minus1);

	/* The published projection must be short. The constant coefficient is a
	 * property of the batch-wide h, so it is checked once, by the caller. */
	result &= proj_norm2_leq(pi->zp, (uint64_t) 4 * PROJ * SIGMA_P * SIGMA_P);

	for (int i = 0; i < SLOTS; i++) {
		inner(u[i], key->b2[i], z, LNP_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(tmp, d[k], com->c2[i][k], k);
			nmod_poly_sub(u[i][k], u[i][k], tmp);
		}
	}
	/* rho * [ sigma(u0) (u0 + d) + sigma(d) u1 ], the part of the quadratic
	 * identity that stays per message. */
	lnp_auto_crt(su, u[SLOT_S], minus1);
	for (int k = 0; k < NCRT; k++) {
		pcrt_poly_mulmod(tmp, d[k], one[k], k);
		nmod_poly_add(lhs[k], u[SLOT_S][k], tmp);
		pcrt_poly_mulmod(lhs[k], su[k], lhs[k], k);
		pcrt_poly_mulmod(tmp, sd[k], u[SLOT_F][k], k);
		nmod_poly_add(lhs[k], lhs[k], tmp);
		pcrt_poly_mulmod(lhs[k], rho[k], lhs[k], k);
		nmod_poly_add(batch->acc[k], batch->acc[k], lhs[k]);
	}

	/* This message's share of sum_l [<B_lj, z_l> - c (T_lj - Z_lj)], which the
	 * batch check compares against v_j once every message has contributed. */
	for (int j = 0; j < LNP_LAMBDA; j++) {
		range_row(b, t, key, com, ctx->P[j], ctx->M[j], ctx->nu[j]);
		inner(lhs, b, z, LNP_WIDTH);
		scalar_crt(zc, ctx->Z[j]);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_sub(t[k], t[k], zc[k]);
			pcrt_poly_mulmod(tmp, d[k], t[k], k);
			nmod_poly_sub(lhs[k], lhs[k], tmp);
			nmod_poly_add(acc[j][k], acc[j][k], lhs[k]);
		}
	}

	nmod_poly_clear(tmp);
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(b[i][k]);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(sd[k]);
		nmod_poly_clear(su[k]);
		nmod_poly_clear(lhs[k]);
		nmod_poly_clear(one[k]);
		nmod_poly_clear(t[k]);
		nmod_poly_clear(zc[k]);
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(u[i][k]);
		}
	}
	return result;
}

void lnp_maskkey_init(lnpmaskkey_t *key) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < MASK_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(key->B1[i][j][k], MODP);
			}
		}
	}
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int j = 0; j < MASK_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(key->b2[i][j][k], MODP);
			}
		}
	}
}

void lnp_maskkey_free(lnpmaskkey_t *key) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < MASK_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(key->B1[i][j][k]);
			}
		}
	}
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int j = 0; j < MASK_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(key->b2[i][j][k]);
			}
		}
	}
}

void lnp_maskkey_gen(lnpmaskkey_t *key, flint_rand_t rand) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < MASK_WIDTH; j++) {
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
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int j = 0; j < MASK_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_zero(key->b2[i][j][k]);
				if (j == HEIGHT + i) {
					nmod_poly_set_coeff_ui(key->b2[i][j][k], 0, 1);
				} else if (j >= HEIGHT + LNP_LAMBDA) {
					commit_sample_rand(key->b2[i][j][k], rand, DEGCRT);
				}
			}
		}
	}
}

void lnp_maskcom_init(lnpmaskcom_t *com) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(com->c1[i][k], MODP);
		}
	}
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(com->c2[i][k], MODP);
		}
	}
}

void lnp_maskcom_free(lnpmaskcom_t *com) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(com->c1[i][k]);
		}
	}
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(com->c2[i][k]);
		}
	}
}

void lnp_mask_commit(lnpmaskcom_t *com, pcrt_poly_t g[LNP_LAMBDA],
		lnpmaskkey_t *key, pcrt_poly_t r[MASK_WIDTH]) {
	for (int i = 0; i < HEIGHT; i++) {
		inner(com->c1[i], key->B1[i], r, MASK_WIDTH);
	}
	for (int i = 0; i < LNP_LAMBDA; i++) {
		inner(com->c2[i], key->b2[i], r, MASK_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_add(com->c2[i][k], com->c2[i][k], g[i][k]);
		}
	}
}

void lnp_garbkey_init(lnpgarbkey_t *key) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < GARB_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(key->B1[i][j][k], MODP);
			}
		}
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < GARB_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(key->b2[i][j][k], MODP);
			}
		}
	}
}

void lnp_garbkey_free(lnpgarbkey_t *key) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < GARB_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(key->B1[i][j][k]);
			}
		}
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < GARB_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(key->b2[i][j][k]);
			}
		}
	}
}

void lnp_garbkey_gen(lnpgarbkey_t *key, flint_rand_t rand) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < GARB_WIDTH; j++) {
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
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < GARB_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_zero(key->b2[i][j][k]);
				if (j == HEIGHT + i) {
					nmod_poly_set_coeff_ui(key->b2[i][j][k], 0, 1);
				} else if (j >= HEIGHT + 2) {
					commit_sample_rand(key->b2[i][j][k], rand, DEGCRT);
				}
			}
		}
	}
}

void lnp_garbcom_init(lnpgarbcom_t *com) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(com->c1[i][k], MODP);
		}
	}
	for (int i = 0; i < 2; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(com->c2[i][k], MODP);
		}
	}
}

void lnp_garbcom_free(lnpgarbcom_t *com) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(com->c1[i][k]);
		}
	}
	for (int i = 0; i < 2; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(com->c2[i][k]);
		}
	}
}

void lnp_garb_commit(lnpgarbcom_t *com, pcrt_poly_t g[2], lnpgarbkey_t *key,
		pcrt_poly_t r[GARB_WIDTH]) {
	for (int i = 0; i < HEIGHT; i++) {
		inner(com->c1[i], key->B1[i], r, GARB_WIDTH);
	}
	for (int i = 0; i < 2; i++) {
		inner(com->c2[i], key->b2[i], r, GARB_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_add(com->c2[i][k], com->c2[i][k], g[i][k]);
		}
	}
}

/* The garbage commitment's Ajtai first message, and the two mask terms of the
 * batched challenge-free term: T gets V_2 + sigma(V_3) once for the batch,
 * where V_i is the mask's contribution to garbage row i. */
void lnp_garb_first(lnpbatch_t *batch, lnpgarbkey_t *gkey,
		pcrt_poly_t yg[GARB_WIDTH]) {
	pcrt_poly_t v2, v3;
	const slong minus1 = 2 * DEGREE - 1;

	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(v2[k], MODP);
		nmod_poly_init(v3[k], MODP);
	}
	for (int i = 0; i < HEIGHT; i++) {
		inner(batch->gw[i], gkey->B1[i], yg, GARB_WIDTH);
	}
	inner(v2, gkey->b2[0], yg, GARB_WIDTH);
	inner(v3, gkey->b2[1], yg, GARB_WIDTH);
	lnp_auto_crt(v3, v3, minus1);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_add(batch->T[k], batch->T[k], v2[k]);
		nmod_poly_add(batch->T[k], batch->T[k], v3[k]);
		nmod_poly_clear(v2[k]);
		nmod_poly_clear(v3[k]);
	}
}

void lnp_batch_init(lnpbatch_t *b) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(b->w[i][k], MODP);
			nmod_poly_init(b->gw[i][k], MODP);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(b->T[k], MODP);
		nmod_poly_init(b->acc[k], MODP);
	}
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(b->h[i][k], MODP);
			nmod_poly_init(b->v[i][k], MODP);
		}
	}
}

void lnp_batch_free(lnpbatch_t *b) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(b->w[i][k]);
			nmod_poly_clear(b->gw[i][k]);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(b->T[k]);
		nmod_poly_clear(b->acc[k]);
	}
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(b->h[i][k]);
			nmod_poly_clear(b->v[i][k]);
		}
	}
}

void lnp_batch_zero(lnpbatch_t *b) {
	for (int i = 0; i < HEIGHT; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_zero(b->w[i][k]);
			nmod_poly_zero(b->gw[i][k]);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_zero(b->T[k]);
		nmod_poly_zero(b->acc[k]);
	}
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_zero(b->h[i][k]);
			nmod_poly_zero(b->v[i][k]);
		}
	}
}

/* The mask commitment's Ajtai first message <B1[i], y_mask>, and its
 * contribution to v, <b2[j], y_mask>. */
void lnp_batch_first(lnpbatch_t *batch, lnpmaskkey_t *mkey,
		pcrt_poly_t ym[MASK_WIDTH]) {
	pcrt_poly_t part;

	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(part[k], MODP);
	}
	for (int i = 0; i < HEIGHT; i++) {
		inner(batch->w[i], mkey->B1[i], ym, MASK_WIDTH);
	}
	for (int j = 0; j < LNP_LAMBDA; j++) {
		inner(part, mkey->b2[j], ym, MASK_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_add(batch->v[j][k], batch->v[j][k], part[k]);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(part[k]);
	}
}

/*
 * The batch check. The first messages w, gw, v and T are not transmitted, so
 * what used to be four comparisons here are now the definitions that recover
 * them; the comparison happens once, when the caller rebuilds the digest over
 * them. What remains a check here is what the digest cannot catch: the norm
 * bounds, and that each aggregated value has zero constant coefficient.
 *
 * Every message has contributed its share of
 * sum_l [<B_lj, z_l> - c (T_lj - Z_lj)] to acc, so what remains is the mask
 * commitment's own term and the comparison against v_j.
 */
int lnp_batch_check(lnpbatch_t *batch, lnpmaskcom_t *mcom, lnpmaskkey_t *mkey,
		lnpgarbcom_t *gcom, lnpgarbkey_t *gkey, pcrt_poly_t d,
		pcrt_poly_t zm[MASK_WIDTH], pcrt_poly_t zg[GARB_WIDTH],
		pcrt_poly_t acc[LNP_LAMBDA]) {
	pcrt_poly_t lhs, t;
	nmod_poly_t tmp, rec;
	int result = 1;

	nmod_poly_init(tmp, MODP);
	nmod_poly_init(rec, MODP);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(lhs[k], MODP);
		nmod_poly_init(t[k], MODP);
	}

	/* The mask opening has to be short and has to open the mask commitment,
	 * or the rows below constrain nothing: they are LNP_LAMBDA equations in
	 * MASK_WIDTH unknowns, so a prover free to choose z_mask could pick any h
	 * with zero constant coefficient and solve for it. Together these two are
	 * the MSIS binding of the mask commitment. */
	for (int i = 0; i < MASK_WIDTH; i++) {
		pcrt_poly_rec(rec, zm[i]);
		result &= commit_norm2_leq(rec,
				(uint64_t) 4 * DEGREE * SIGMA_B * SIGMA_B);
	}
	for (int i = 0; i < HEIGHT; i++) {
		inner(lhs, mkey->B1[i], zm, MASK_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(tmp, d[k], mcom->c1[i][k], k);
			nmod_poly_sub(batch->w[i][k], lhs[k], tmp);
		}
	}

	/* The garbage commitment is opened the same way, and for the same reason:
	 * the batched quadratic identity below is one equation in GARB_WIDTH
	 * unknowns without it. */
	for (int i = 0; i < GARB_WIDTH; i++) {
		pcrt_poly_rec(rec, zg[i]);
		result &= commit_norm2_leq(rec,
				(uint64_t) 4 * DEGREE * SIGMA_B * SIGMA_B);
	}
	for (int i = 0; i < HEIGHT; i++) {
		inner(lhs, gkey->B1[i], zg, GARB_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(tmp, d[k], gcom->c1[i][k], k);
			nmod_poly_sub(batch->gw[i][k], lhs[k], tmp);
		}
	}

	/* The batched quadratic relation: every message has added its rho-weighted
	 * share to acc, and what remains is the garbage commitment's own two rows,
	 * U2 + sigma(U3), against the batched challenge-free term T. */
	{
		pcrt_poly_t u2, u3;
		const slong minus1 = 2 * DEGREE - 1;

		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(u2[k], MODP);
			nmod_poly_init(u3[k], MODP);
		}
		inner(u2, gkey->b2[0], zg, GARB_WIDTH);
		inner(u3, gkey->b2[1], zg, GARB_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(tmp, d[k], gcom->c2[0][k], k);
			nmod_poly_sub(u2[k], u2[k], tmp);
			pcrt_poly_mulmod(tmp, d[k], gcom->c2[1][k], k);
			nmod_poly_sub(u3[k], u3[k], tmp);
		}
		lnp_auto_crt(u3, u3, minus1);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_add(lhs[k], batch->acc[k], u2[k]);
			nmod_poly_add(lhs[k], lhs[k], u3[k]);
			nmod_poly_set(batch->T[k], lhs[k]);
			nmod_poly_clear(u2[k]);
			nmod_poly_clear(u3[k]);
		}
	}

	/* The constant coefficient of each aggregated value, which is the whole
	 * point of the construction. */
	for (int j = 0; j < LNP_LAMBDA; j++) {
		pcrt_poly_rec(rec, batch->h[j]);
		result &= (nmod_poly_get_coeff_ui(rec, 0) == 0);
	}

	for (int j = 0; j < LNP_LAMBDA; j++) {
		inner(lhs, mkey->b2[j], zm, MASK_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_sub(t[k], mcom->c2[j][k], batch->h[j][k]);
			pcrt_poly_mulmod(tmp, d[k], t[k], k);
			nmod_poly_sub(lhs[k], lhs[k], tmp);
			nmod_poly_add(lhs[k], lhs[k], acc[j][k]);
			nmod_poly_set(batch->v[j][k], lhs[k]);
		}
	}

	nmod_poly_clear(tmp);
	nmod_poly_clear(rec);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(lhs[k]);
		nmod_poly_clear(t[k]);
	}
	return result;
}

/*============================================================================*/
/* Tests                                                                      */
/*============================================================================*/

#ifdef LNP_MAIN

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















/* Reconstruct the constant coefficient of a CRT-represented polynomial. */
static ulong const_coeff(pcrt_poly_t a) {
	nmod_poly_t t;
	ulong r;

	nmod_poly_init(t, MODP);
	pcrt_poly_rec(t, a);
	r = nmod_poly_get_coeff_ui(t, 0);
	nmod_poly_clear(t);
	return r;
}



/* Drive the batched proof the way the shuffle does, with a batch of one
 * message: commit the masks, set up, first message, one challenge, one
 * response, then the per-message and batch checks. Tests go through this so
 * that they exercise the same path as the real caller. */
static lnpmaskkey_t tst_mkey;
static lnpmaskcom_t tst_mcom;
static lnpgarbkey_t tst_gkey;
static lnpgarbcom_t tst_gcom;
static pcrt_poly_t tst_zg[GARB_WIDTH], tst_rg[GARB_WIDTH], tst_rho;
static lnpbatch_t tst_batch;
static pcrt_poly_t tst_z[LNP_WIDTH], tst_zm[MASK_WIDTH], tst_ajtai[HEIGHT];
static uint8_t tst_digest[SHA256HashSize];
static pcrt_poly_t tst_mr[MASK_WIDTH];
static int tst_ready = 0;

static void tst_setup(flint_rand_t rng) {
	if (tst_ready) {
		return;
	}
	lnp_maskkey_init(&tst_mkey);
	lnp_maskkey_gen(&tst_mkey, rng);
	lnp_maskcom_init(&tst_mcom);
	lnp_garbkey_init(&tst_gkey);
	lnp_garbkey_gen(&tst_gkey, rng);
	lnp_garbcom_init(&tst_gcom);
	for (int i = 0; i < GARB_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(tst_zg[i][k], MODP);
			nmod_poly_init(tst_rg[i][k], MODP);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(tst_rho[k], MODP);
	}
	lnp_batch_init(&tst_batch);
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(tst_z[i][k], MODP);
		}
	}
	for (int i = 0; i < MASK_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(tst_zm[i][k], MODP);
			nmod_poly_init(tst_mr[i][k], MODP);
		}
	}
	for (int i = 0; i < HEIGHT; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(tst_ajtai[i][k], MODP);
		}
	}
	tst_ready = 1;
}

/* The challenge as a function of the digest, so that the harness can derive it
 * before the first messages exist, exactly as the shuffle's verifier does. */
static void challenge_from_digest(pcrt_poly_t d,
		const uint8_t hash[SHA256HashSize]) {
	nmod_poly_t c;
	uint32_t buf;

	nmod_poly_init(c, MODP);
	fastrandombytes_setseed((uint8_t *) hash);
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

static void bin_local_hash(pcrt_poly_t d, lnpkey_t *key, lnpcom_t *com,
		lnpbinproof_t *pi, pcrt_poly_t ajtai[HEIGHT], uint8_t *digest) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	uint32_t buf;
	nmod_poly_t c;

	SHA256Reset(&sha);
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
			hash_poly(&sha, ajtai[i][k]);
			hash_poly(&sha, tst_mcom.c1[i][k]);
			hash_poly(&sha, tst_batch.w[i][k]);
		}
		for (int i = 0; i < SLOTS; i++) {
			hash_poly(&sha, com->c2[i][k]);
		}
		for (int i = 0; i < LNP_LAMBDA; i++) {
			hash_poly(&sha, tst_mcom.c2[i][k]);
			hash_poly(&sha, tst_batch.h[i][k]);
			hash_poly(&sha, tst_batch.v[i][k]);
		}
		for (int i = 0; i < HEIGHT; i++) {
			hash_poly(&sha, tst_gcom.c1[i][k]);
			hash_poly(&sha, tst_batch.gw[i][k]);
		}
		for (int i = 0; i < 2; i++) {
			hash_poly(&sha, tst_gcom.c2[i][k]);
		}
		hash_poly(&sha, tst_batch.T[k]);
	}
	SHA256Input(&sha, (const uint8_t *)pi->zp, PROJ * sizeof(ulong));
	SHA256Result(&sha, hash);
	if (digest != NULL) {
		memcpy(digest, hash, SHA256HashSize);
	}

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

/* The batching challenge. It has to be drawn after every message's commitment,
 * because the garbage terms are weighted by it, and before the garbage terms
 * are committed, because a prover that knew it could adapt them. */
static void tst_rho_hash(pcrt_poly_t rho, lnpkey_t *key, lnpcom_t *com) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	nmod_poly_t c;

	SHA256Reset(&sha);
	for (int k = 0; k < NCRT; k++) {
		for (int i = 0; i < HEIGHT; i++) {
			hash_poly(&sha, com->c1[i][k]);
			hash_poly(&sha, tst_mcom.c1[i][k]);
		}
		for (int i = 0; i < SLOTS; i++) {
			hash_poly(&sha, com->c2[i][k]);
		}
		for (int i = 0; i < LNP_LAMBDA; i++) {
			hash_poly(&sha, tst_mcom.c2[i][k]);
		}
		for (int i = 0; i < LNP_WIDTH; i++) {
			hash_poly(&sha, key->b2[0][i][k]);
		}
	}
	SHA256Result(&sha, hash);

	/* Same shape as the opening challenge: NONZERO ones in DEGREE positions.
	 * Differences of distinct challenges are then short enough for Lemma 1 to
	 * make them invertible, which is what the aggregation needs. */
	nmod_poly_init(c, MODP);
	fastrandombytes_setseed(hash);
	nmod_poly_fit_length(c, DEGREE);
	for (int i = 0; i < NONZERO; i++) {
		uint32_t buf;

		fastrandombytes((unsigned char *)&buf, sizeof(buf));
		buf = buf % DEGREE;
		while (nmod_poly_get_coeff_ui(c, buf) != 0) {
			fastrandombytes((unsigned char *)&buf, sizeof(buf));
			buf = buf % DEGREE;
		}
		nmod_poly_set_coeff_ui(c, buf, 1);
	}
	pcrt_poly_reduce(rho[0], c, 0);
	pcrt_poly_reduce(rho[1], c, 1);
	nmod_poly_clear(c);
}

static int bin_prove_local(lnpbinproof_t *pi, lnpbinctx_t *ctx, lnpcom_t *com,
		pcrt_poly_t s, pcrt_poly_t f, nmod_poly_t w, pcrt_poly_t g[LNP_LAMBDA],
		lnpkey_t *key, pcrt_poly_t r[LNP_WIDTH], flint_rand_t rng) {
	pcrt_poly_t y[LNP_WIDTH], cr[LNP_WIDTH], ym[MASK_WIDTH], cm[MASK_WIDTH], d;
	pcrt_poly_t yg[GARB_WIDTH], cg[GARB_WIDTH], garb[2];
	int64_t dot, norm;
	int rej, ok;

	tst_setup(rng);
	for (int i = 0; i < MASK_WIDTH; i++) {
		commit_sample_short_crt(tst_mr[i]);
	}
	for (int i = 0; i < GARB_WIDTH; i++) {
		commit_sample_short_crt(tst_rg[i]);
	}
	lnp_mask_commit(&tst_mcom, g, &tst_mkey, tst_mr);

	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(y[i][k], MODP);
			nmod_poly_init(cr[i][k], MODP);
		}
	}
	for (int i = 0; i < MASK_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(ym[i][k], MODP);
			nmod_poly_init(cm[i][k], MODP);
		}
	}
	for (int i = 0; i < GARB_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(yg[i][k], MODP);
			nmod_poly_init(cg[i][k], MODP);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(d[k], MODP);
		nmod_poly_init(garb[0][k], MODP);
		nmod_poly_init(garb[1][k], MODP);
	}

	lnp_batch_zero(&tst_batch);
	ok = lnp_bin_setup(pi, ctx, &tst_batch, com, &tst_mcom, s, f, w, key);
	/* The shared mask is added once, not once per message. */
	for (int j = 0; j < LNP_LAMBDA; j++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_add(tst_batch.h[j][k], tst_batch.h[j][k], g[j][k]);
		}
	}

	/* The batching challenge comes after every message's commitment and before
	 * the garbage terms, which depend on it. With a batch of one there is
	 * still a rho, so that this path exercises what the shuffle does. */
	tst_rho_hash(tst_rho, key, com);

	do {
		for (int j = 0; j < LNP_LAMBDA; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_zero(tst_batch.v[j][k]);
			}
		}
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_zero(tst_batch.T[k]);
			nmod_poly_zero(garb[0][k]);
			nmod_poly_zero(garb[1][k]);
		}
		for (int i = 0; i < LNP_WIDTH; i++) {
			commit_sample_gauss_batch_crt(y[i]);
		}
		for (int i = 0; i < MASK_WIDTH; i++) {
			commit_sample_gauss_batch_crt(ym[i]);
		}
		for (int i = 0; i < GARB_WIDTH; i++) {
			commit_sample_gauss_batch_crt(yg[i]);
		}
		lnp_bin_first(ctx, &tst_batch, com, s, key, r, y, tst_rho, garb);
		lnp_garb_commit(&tst_gcom, garb, &tst_gkey, tst_rg);
		lnp_garb_first(&tst_batch, &tst_gkey, yg);
		lnp_batch_first(&tst_batch, &tst_mkey, ym);
		for (int i = 0; i < HEIGHT; i++) {
			inner(tst_ajtai[i], key->B1[i], y, LNP_WIDTH);
		}
		bin_local_hash(d, key, com, pi, tst_ajtai, tst_digest);
		dot = norm = 0;
		for (int i = 0; i < LNP_WIDTH; i++) {
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_mulmod(cr[i][k], d[k], r[i][k], k);
				nmod_poly_add(tst_z[i][k], y[i][k], cr[i][k]);
			}
		}
		for (int i = 0; i < MASK_WIDTH; i++) {
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_mulmod(cm[i][k], d[k], tst_mr[i][k], k);
				nmod_poly_add(tst_zm[i][k], ym[i][k], cm[i][k]);
			}
		}
		for (int i = 0; i < GARB_WIDTH; i++) {
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_mulmod(cg[i][k], d[k], tst_rg[i][k], k);
				nmod_poly_add(tst_zg[i][k], yg[i][k], cg[i][k]);
			}
		}
		commit_rej_accumulate(&dot, &norm, tst_z, cr, LNP_WIDTH);
		commit_rej_accumulate(&dot, &norm, tst_zm, cm, MASK_WIDTH);
		commit_rej_accumulate(&dot, &norm, tst_zg, cg, GARB_WIDTH);
		rej = commit_rej_decide(dot, norm, (uint64_t) SIGMA_B * SIGMA_B);
	} while (rej);

	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(y[i][k]);
			nmod_poly_clear(cr[i][k]);
		}
	}
	for (int i = 0; i < MASK_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(ym[i][k]);
			nmod_poly_clear(cm[i][k]);
		}
	}
	for (int i = 0; i < GARB_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(yg[i][k]);
			nmod_poly_clear(cg[i][k]);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(d[k]);
		nmod_poly_clear(garb[0][k]);
		nmod_poly_clear(garb[1][k]);
	}
	return ok;
}

static int bin_verify_local(lnpbinproof_t *pi, lnpbinctx_t *ctx, lnpcom_t *com,
		lnpkey_t *key) {
	pcrt_poly_t d, acc[LNP_LAMBDA], lhs, rhs;
	nmod_poly_t rec, tmp;
	int result = 1;

	nmod_poly_init(rec, MODP);
	nmod_poly_init(tmp, MODP);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(d[k], MODP);
		nmod_poly_init(lhs[k], MODP);
		nmod_poly_init(rhs[k], MODP);
	}
	for (int j = 0; j < LNP_LAMBDA; j++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(acc[j][k], MODP);
		}
	}
	lnp_bin_public(ctx, pi, com, &tst_mcom, key);
	/* The first messages are no longer compared inside lnp_batch_check; that
	 * function rebuilds them. So this has to close the loop the same way
	 * shuffle_verifier does, or the checks it replaced are simply gone. The
	 * challenge comes from the prover's digest, the checks below rebuild
	 * w, gw, v and T, and the digest is rebuilt over them at the end. */
	challenge_from_digest(d, tst_digest);
	for (int i = 0; i < LNP_WIDTH; i++) {
		pcrt_poly_rec(rec, tst_z[i]);
		result &= commit_norm2_leq(rec,
				(uint64_t) 4 * DEGREE * SIGMA_B * SIGMA_B);
	}
	for (int i = 0; i < HEIGHT; i++) {
		inner(lhs, key->B1[i], tst_z, LNP_WIDTH);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(tmp, d[k], com->c1[i][k], k);
			nmod_poly_add(rhs[k], tst_ajtai[i][k], tmp);
			result &= nmod_poly_equal(lhs[k], rhs[k]);
		}
	}
	/* The verifier zeroes the quadratic accumulator, then every message adds
	 * its rho-weighted share before the batch check settles it. */
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_zero(tst_batch.acc[k]);
	}
	result &= lnp_bin_check(pi, ctx, &tst_batch, com, key, d, tst_z, acc,
			tst_rho);
	result &= lnp_batch_check(&tst_batch, &tst_mcom, &tst_mkey, &tst_gcom,
			&tst_gkey, d, tst_zm, tst_zg, acc);
	{
		uint8_t rebuilt[SHA256HashSize];
		pcrt_poly_t ignored;

		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(ignored[k], MODP);
		}
		bin_local_hash(ignored, key, com, pi, tst_ajtai, rebuilt);
		result &= (memcmp(rebuilt, tst_digest, SHA256HashSize) == 0);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(ignored[k]);
		}
	}

	nmod_poly_clear(rec);
	nmod_poly_clear(tmp);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(d[k]);
		nmod_poly_clear(lhs[k]);
		nmod_poly_clear(rhs[k]);
	}
	for (int j = 0; j < LNP_LAMBDA; j++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(acc[j][k]);
		}
	}
	return result;
}

/* The two halves composed: the product relation says slot 1 holds
 * sigma_{-1}(s) (s - ones), and the constant-coefficient proof says its
 * constant coefficient is zero. Together they say that every coefficient of s
 * is 0 or 1, given the norm bound that B5 would supply. */
static int isbin_full(lnpkey_t *key, flint_rand_t rng, nmod_poly_t s) {
	lnpcom_t com;
	lnpproof_t pi;
	lnpbinproof_t pc;
	lnpbinctx_t ctx_pc;
	nmod_poly_t wraw;
	pcrt_poly_t wc;
	pcrt_poly_t r[LNP_WIDTH], m[SLOTS], g[LNP_LAMBDA];
	int result;

	lnp_com_init(&com);
	lnp_proof_init(&pi);
	lnp_binproof_init(&pc);
	lnp_binctx_init(&ctx_pc);
	nmod_poly_init(wraw, MODP);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(wc[k], MODP);
	}
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(r[i][k], MODP);
		}
		commit_sample_short_crt(r[i]);
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(m[i][k], MODP);
			nmod_poly_zero(m[i][k]);
		}
	}
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(g[i][k], MODP);
		}
		lnp_sample_ct_zero(g[i], rng);
	}
	pcrt_poly_reduce(m[0][0], s, 0);
	pcrt_poly_reduce(m[0][1], s, 1);
	lnp_isbin_product(m[SLOT_F], m[0]);
	lnp_commit(&com, m, key, r);

	/* The projection can be rejected, and the honest prover recommits and
	 * tries again rather than sending the rejected attempt. */
	for (int tries = 0; tries < 64; tries++) {
		lnp_sample_proj_mask(wc, wraw);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_set(m[SLOT_W][k], wc[k]);
		}
		lnp_commit(&com, m, key, r);
		if (bin_prove_local(&pc, &ctx_pc, &com, m[0], m[SLOT_F], wraw, g, key,
				r, rng)) {
			break;
		}
	}
	result = bin_verify_local(&pc, &ctx_pc, &com, key);

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
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(g[i][k]);
		}
	}
	lnp_com_free(&com);
	lnp_proof_free(&pi);
	lnp_binproof_free(&pc);
	lnp_binctx_free(&ctx_pc);
	nmod_poly_clear(wraw);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(wc[k]);
	}
	return result;
}



/* Select which phases to run: "test", "bench", or neither for both. Same
 * helper as in commit.c and shuffle.c; it stays local to each binary so that
 * the test harness does not have to be shared between them. */
static int phase_selected(int argc, char *argv[], const char *phase) {
	return argc < 2 || strcmp(argv[1], phase) == 0;
}


/* A uniformly random binary witness, in both representations. */
static void binary_witness(pcrt_poly_t out, nmod_poly_t raw) {
	uint64_t buf = 0;

	nmod_poly_zero(raw);
	nmod_poly_fit_length(raw, DEGREE);
	for (int i = 0; i < DEGREE; i++) {
		if (i % 64 == 0) {
			getrandom(&buf, sizeof(buf), 0);
		}
		nmod_poly_set_coeff_ui(raw, i, (buf >> (i % 64)) & 1);
	}
	pcrt_poly_reduce(out[0], raw, 0);
	pcrt_poly_reduce(out[1], raw, 1);
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


static void test_isbin(flint_rand_t rng) {
	lnpkey_t key;
	lnpcom_t com;
	lnpbinproof_t pi;
	lnpbinctx_t ctx_pi;
	pcrt_poly_t r[LNP_WIDTH], m[SLOTS], g[LNP_LAMBDA], wc;
	nmod_poly_t wraw;
	nmod_poly_t s;
	uint64_t buf;

	lnp_keyinit(&key);
	lnp_com_init(&com);
	lnp_binproof_init(&pi);
	lnp_binctx_init(&ctx_pi);
	nmod_poly_init(wraw, MODP);
	lnp_keygen(&key, rng);
	nmod_poly_init(s, MODP);
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(r[i][k], MODP);
		}
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(m[i][k], MODP);
			nmod_poly_zero(m[i][k]);
		}
	}
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(g[i][k], MODP);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(wc[k], MODP);
	}

	TEST_BEGIN("constant coefficient of the product is zero for binary s") {
		nmod_poly_zero(s);
		nmod_poly_fit_length(s, DEGREE);
		for (int i = 0; i < DEGREE; i++) {
			if (i % 64 == 0) {
				getrandom(&buf, sizeof(buf), 0);
			}
			nmod_poly_set_coeff_ui(s, i, (buf >> (i % 64)) & 1);
		}
		pcrt_poly_reduce(m[0][0], s, 0);
		pcrt_poly_reduce(m[0][1], s, 1);
		lnp_isbin_product(m[1], m[0]);
		TEST_ASSERT(const_coeff(m[1]) == 0, end);
	} TEST_END;

	TEST_BEGIN("constant coefficient is non-zero as soon as s is not binary") {
		nmod_poly_zero(s);
		nmod_poly_fit_length(s, DEGREE);
		for (int i = 0; i < DEGREE; i++) {
			if (i % 64 == 0) {
				getrandom(&buf, sizeof(buf), 0);
			}
			nmod_poly_set_coeff_ui(s, i, (buf >> (i % 64)) & 1);
		}
		/* A single coefficient of 2 contributes 2 * 1 = 2 to the sum, far
		 * below the modulus, so the sum cannot wrap back to zero. */
		nmod_poly_set_coeff_ui(s, 7, 2);
		pcrt_poly_reduce(m[0][0], s, 0);
		pcrt_poly_reduce(m[0][1], s, 1);
		lnp_isbin_product(m[1], m[0]);
		TEST_ASSERT(const_coeff(m[1]) != 0, end);
	} TEST_END;

	TEST_BEGIN("is_bin product proof accepts an honest product") {
		for (int i = 0; i < LNP_WIDTH; i++) {
			commit_sample_short_crt(r[i]);
		}
		binary_witness(m[SLOT_S], s);
		lnp_isbin_product(m[SLOT_F], m[SLOT_S]);
		for (int i = 0; i < LNP_LAMBDA; i++) {
			lnp_sample_ct_zero(g[i], rng);
		}
		/* The projection can be rejected, and the honest prover recommits and
		 * tries again rather than sending the rejected attempt. */
		int tries;

		for (tries = 0; tries < 64; tries++) {
			lnp_sample_proj_mask(wc, wraw);
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_set(m[SLOT_W][k], wc[k]);
			}
			lnp_commit(&com, m, &key, r);
			if (bin_prove_local(&pi, &ctx_pi, &com, m[SLOT_S], m[SLOT_F],
					wraw, g, &key, r, rng)) {
				break;
			}
		}
		TEST_ASSERT(tries < 64, end);
		TEST_ASSERT(bin_verify_local(&pi, &ctx_pi, &com, &key) == 1, end);
	} TEST_END;

	TEST_BEGIN("is_bin product proof rejects a wrong product") {
		for (int i = 0; i < LNP_WIDTH; i++) {
			commit_sample_short_crt(r[i]);
		}
		/* A binary witness, so the constant coefficient and the norm bound
		 * are both satisfied, but a product that is off by one. Only the
		 * product relation can reject this. */
		binary_witness(m[SLOT_S], s);
		lnp_isbin_product(m[SLOT_F], m[SLOT_S]);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_set_coeff_ui(m[SLOT_F][k], 1,
					nmod_add(nmod_poly_get_coeff_ui(m[SLOT_F][k], 1), 1,
					m[SLOT_F][k]->mod));
		}
		for (int i = 0; i < LNP_LAMBDA; i++) {
			lnp_sample_ct_zero(g[i], rng);
		}
		/* The projection can be rejected, and the honest prover recommits and
		 * tries again rather than sending the rejected attempt. */
		int tries;

		for (tries = 0; tries < 64; tries++) {
			lnp_sample_proj_mask(wc, wraw);
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_set(m[SLOT_W][k], wc[k]);
			}
			lnp_commit(&com, m, &key, r);
			if (bin_prove_local(&pi, &ctx_pi, &com, m[SLOT_S], m[SLOT_F],
					wraw, g, &key, r, rng)) {
				break;
			}
		}
		TEST_ASSERT(tries < 64, end);
		TEST_ASSERT(bin_verify_local(&pi, &ctx_pi, &com, &key) == 0, end);
	} TEST_END;

  end:
	nmod_poly_clear(s);
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
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(g[i][k]);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(wc[k]);
	}
	lnp_keyfree(&key);
	lnp_com_free(&com);
	lnp_binproof_free(&pi);
	lnp_binctx_free(&ctx_pi);
	nmod_poly_clear(wraw);
}

static void test_adaptive(flint_rand_t rng) {
	lnpkey_t key;
	lnpcom_t com;
	lnpbinproof_t pi;
	lnpbinctx_t ctx_pi;
	nmod_poly_t wraw;
	pcrt_poly_t wc;
	pcrt_poly_t r[LNP_WIDTH], m[SLOTS], g[LNP_LAMBDA];

	lnp_keyinit(&key);
	lnp_com_init(&com);
	lnp_binproof_init(&pi);
	lnp_binctx_init(&ctx_pi);
	nmod_poly_init(wraw, MODP);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(wc[k], MODP);
	}
	lnp_keygen(&key, rng);
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
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(g[i][k], MODP);
		}
	}

	TEST_BEGIN("a mask chosen to cancel the constant term does not help") {
		nmod_poly_t tt, tw;
		ulong nu_seen[LNP_LAMBDA];

		nmod_poly_init(tt, MODP);
		nmod_poly_init(tw, MODP);
		for (int i = 0; i < LNP_WIDTH; i++) {
			commit_sample_short_crt(r[i]);
		}
		for (int i = 0; i < SLOTS; i++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_zero(m[i][k]);
			}
		}
		/* A witness with a single coefficient of 2, so the product's constant
		 * coefficient is 2 * 1 = 2: the product relation holds and only the
		 * constant-coefficient check can reject. */
		binary_witness(m[SLOT_S], tw);
		nmod_poly_set_coeff_ui(tw, 9, 2);
		pcrt_poly_reduce(m[SLOT_S][0], tw, 0);
		pcrt_poly_reduce(m[SLOT_S][1], tw, 1);
		lnp_isbin_product(m[SLOT_F], m[SLOT_S]);
		for (int i = 0; i < LNP_LAMBDA; i++) {
			lnp_sample_ct_zero(g[i], rng);
		}
		/* The projection can be rejected, and the honest prover recommits and
		 * tries again rather than sending the rejected attempt. */
		int tries;

		for (tries = 0; tries < 64; tries++) {
			lnp_sample_proj_mask(wc, wraw);
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_set(m[SLOT_W][k], wc[k]);
			}
			lnp_commit(&com, m, &key, r);
			if (bin_prove_local(&pi, &ctx_pi, &com, m[SLOT_S], m[SLOT_F],
					wraw, g, &key, r, rng)) {
				break;
			}
		}
		TEST_ASSERT(tries < 64, end);
		/* The cheating prover reads the scalars it would face and rewrites
		 * each mask so that nu_j * ct(f) is cancelled. */
		lnp_scalars_for_test(nu_seen, &key, &com, &tst_mcom, pi.zp);
		for (int i = 0; i < LNP_LAMBDA; i++) {
			pcrt_poly_rec(tt, g[i]);
			nmod_poly_set_coeff_ui(tt, 0,
					nmod_neg(nmod_mul(nu_seen[i], 2, tt->mod), tt->mod));
			pcrt_poly_reduce(g[i][0], tt, 0);
			pcrt_poly_reduce(g[i][1], tt, 1);
		}
		/* Committing the adapted masks changes the commitment, so the scalars
		 * it now faces are not the ones it prepared for. */
		lnp_commit(&com, m, &key, r);
		bin_prove_local(&pi, &ctx_pi, &com, m[SLOT_S], m[SLOT_F], wraw, g, &key, r, rng);
		nmod_poly_clear(tt);
		nmod_poly_clear(tw);
		TEST_ASSERT(bin_verify_local(&pi, &ctx_pi, &com, &key) == 0, end);
	} TEST_END;


  end:
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
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(g[i][k]);
		}
	}
	lnp_keyfree(&key);
	lnp_com_free(&com);
	lnp_binproof_free(&pi);
	lnp_binctx_free(&ctx_pi);
	nmod_poly_clear(wraw);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(wc[k]);
	}
}

static void test_isbin_full(flint_rand_t rng) {
	lnpkey_t key;
	nmod_poly_t s;
	uint64_t buf;

	lnp_keyinit(&key);
	lnp_keygen(&key, rng);
	nmod_poly_init(s, MODP);

	TEST_BEGIN("is_bin accepts a binary witness") {
		nmod_poly_zero(s);
		nmod_poly_fit_length(s, DEGREE);
		for (int i = 0; i < DEGREE; i++) {
			if (i % 64 == 0) {
				getrandom(&buf, sizeof(buf), 0);
			}
			nmod_poly_set_coeff_ui(s, i, (buf >> (i % 64)) & 1);
		}
		TEST_ASSERT(isbin_full(&key, rng, s) == 1, end);
	} TEST_END;

	TEST_BEGIN("is_bin rejects a witness with a coefficient of 2") {
		nmod_poly_zero(s);
		nmod_poly_fit_length(s, DEGREE);
		for (int i = 0; i < DEGREE; i++) {
			if (i % 64 == 0) {
				getrandom(&buf, sizeof(buf), 0);
			}
			nmod_poly_set_coeff_ui(s, i, (buf >> (i % 64)) & 1);
		}
		nmod_poly_set_coeff_ui(s, 11, 2);
		TEST_ASSERT(isbin_full(&key, rng, s) == 0, end);
	} TEST_END;

	TEST_BEGIN("is_bin rejects a witness with a coefficient of -1") {
		nmod_poly_zero(s);
		nmod_poly_fit_length(s, DEGREE);
		for (int i = 0; i < DEGREE; i++) {
			if (i % 64 == 0) {
				getrandom(&buf, sizeof(buf), 0);
			}
			nmod_poly_set_coeff_ui(s, i, (buf >> (i % 64)) & 1);
		}
		nmod_poly_set_coeff_ui(s, 23, MODP - 1);
		TEST_ASSERT(isbin_full(&key, rng, s) == 0, end);
	} TEST_END;

  end:
	nmod_poly_clear(s);
	lnp_keyfree(&key);
}

static void test_range(flint_rand_t rng) {
	lnpkey_t key;
	lnpcom_t com;
	lnpbinproof_t pi;
	lnpbinctx_t ctx_pi;
	pcrt_poly_t r[LNP_WIDTH], m[SLOTS], g[LNP_LAMBDA], w;
	nmod_poly_t s, wraw;
	uint64_t buf;
	int tries;

	lnp_keyinit(&key);
	lnp_com_init(&com);
	lnp_binproof_init(&pi);
	lnp_binctx_init(&ctx_pi);
	lnp_keygen(&key, rng);
	nmod_poly_init(s, MODP);
	nmod_poly_init(wraw, MODP);
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
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(g[i][k], MODP);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(w[k], MODP);
	}

	TEST_BEGIN("range proof accepts a binary witness") {
		nmod_poly_zero(s);
		nmod_poly_fit_length(s, DEGREE);
		for (int i = 0; i < DEGREE; i++) {
			if (i % 64 == 0) {
				getrandom(&buf, sizeof(buf), 0);
			}
			nmod_poly_set_coeff_ui(s, i, (buf >> (i % 64)) & 1);
		}
		/* Rejection sampling needs a fresh mask, and the mask is committed,
		 * so a rejected transcript is retried by recommitting. */
		for (tries = 0; tries < 64; tries++) {
			for (int i = 0; i < LNP_WIDTH; i++) {
				commit_sample_short_crt(r[i]);
			}
			for (int i = 0; i < SLOTS; i++) {
				for (int k = 0; k < NCRT; k++) {
					nmod_poly_zero(m[i][k]);
				}
			}
			pcrt_poly_reduce(m[SLOT_S][0], s, 0);
			pcrt_poly_reduce(m[SLOT_S][1], s, 1);
			/* The merged proof also checks the product relation, so slot
			 * SLOT_F has to hold the honest product of the witness. */
			lnp_isbin_product(m[SLOT_F], m[SLOT_S]);
			lnp_sample_proj_mask(w, wraw);
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_set(m[SLOT_W][k], w[k]);
			}
			for (int i = 0; i < LNP_LAMBDA; i++) {
				lnp_sample_ct_zero(g[i], rng);
			}
			lnp_commit(&com, m, &key, r);
			if (bin_prove_local(&pi, &ctx_pi, &com, m[SLOT_S], m[SLOT_F],
					wraw, g, &key, r, rng)) {
				break;
			}
		}
		TEST_ASSERT(tries < 64, end);
		TEST_ASSERT(bin_verify_local(&pi, &ctx_pi, &com, &key) == 1, end);
	} TEST_END;

	TEST_ONCE("range proof rejects a witness that is far too large") {
		/* One coefficient near the modulus makes the projection huge, so the
		 * published vector cannot pass the norm check however it is masked. */
		nmod_poly_zero(s);
		nmod_poly_fit_length(s, DEGREE);
		nmod_poly_set_coeff_ui(s, 3, MODP / 4);
		for (int i = 0; i < LNP_WIDTH; i++) {
			commit_sample_short_crt(r[i]);
		}
		for (int i = 0; i < SLOTS; i++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_zero(m[i][k]);
			}
		}
		pcrt_poly_reduce(m[SLOT_S][0], s, 0);
		pcrt_poly_reduce(m[SLOT_S][1], s, 1);
		lnp_isbin_product(m[SLOT_F], m[SLOT_S]);
		lnp_sample_proj_mask(w, wraw);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_set(m[SLOT_W][k], w[k]);
		}
		for (int i = 0; i < LNP_LAMBDA; i++) {
			lnp_sample_ct_zero(g[i], rng);
		}
		lnp_commit(&com, m, &key, r);
		bin_prove_local(&pi, &ctx_pi, &com, m[SLOT_S], m[SLOT_F], wraw, g, &key, r, rng);
		TEST_ASSERT(bin_verify_local(&pi, &ctx_pi, &com, &key) == 0, end);
	} TEST_END;

	TEST_ONCE("range proof rejects a tampered projection") {
		nmod_poly_zero(s);
		nmod_poly_fit_length(s, DEGREE);
		for (int i = 0; i < DEGREE; i++) {
			nmod_poly_set_coeff_ui(s, i, i & 1);
		}
		for (tries = 0; tries < 64; tries++) {
			for (int i = 0; i < LNP_WIDTH; i++) {
				commit_sample_short_crt(r[i]);
			}
			for (int i = 0; i < SLOTS; i++) {
				for (int k = 0; k < NCRT; k++) {
					nmod_poly_zero(m[i][k]);
				}
			}
			pcrt_poly_reduce(m[SLOT_S][0], s, 0);
			pcrt_poly_reduce(m[SLOT_S][1], s, 1);
			/* The merged proof also checks the product relation, so slot
			 * SLOT_F has to hold the honest product of the witness. */
			lnp_isbin_product(m[SLOT_F], m[SLOT_S]);
			lnp_sample_proj_mask(w, wraw);
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_set(m[SLOT_W][k], w[k]);
			}
			for (int i = 0; i < LNP_LAMBDA; i++) {
				lnp_sample_ct_zero(g[i], rng);
			}
			lnp_commit(&com, m, &key, r);
			if (bin_prove_local(&pi, &ctx_pi, &com, m[SLOT_S], m[SLOT_F],
					wraw, g, &key, r, rng)) {
				break;
			}
		}
		pi.zp[5] = nmod_add(pi.zp[5], 1, s->mod);
		TEST_ASSERT(bin_verify_local(&pi, &ctx_pi, &com, &key) == 0, end);
	} TEST_END;

	/* The mask commitment is opened by z_mask, and both halves of that opening
	 * have to be checked. Without them the batch rows are LNP_LAMBDA equations
	 * in MASK_WIDTH unknowns, and a prover could choose any aggregated values
	 * with zero constant coefficient and solve for z_mask, which would make
	 * the constant-coefficient argument vacuous. */
	TEST_ONCE("is_bin rejects a mask opening that misses the commitment") {
		nmod_poly_zero(s);
		nmod_poly_fit_length(s, DEGREE);
		for (int i = 0; i < DEGREE; i++) {
			nmod_poly_set_coeff_ui(s, i, i & 1);
		}
		for (tries = 0; tries < 64; tries++) {
			for (int i = 0; i < LNP_WIDTH; i++) {
				commit_sample_short_crt(r[i]);
			}
			for (int i = 0; i < SLOTS; i++) {
				for (int k = 0; k < NCRT; k++) {
					nmod_poly_zero(m[i][k]);
				}
			}
			pcrt_poly_reduce(m[SLOT_S][0], s, 0);
			pcrt_poly_reduce(m[SLOT_S][1], s, 1);
			lnp_isbin_product(m[SLOT_F], m[SLOT_S]);
			lnp_sample_proj_mask(w, wraw);
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_set(m[SLOT_W][k], w[k]);
			}
			for (int i = 0; i < LNP_LAMBDA; i++) {
				lnp_sample_ct_zero(g[i], rng);
			}
			lnp_commit(&com, m, &key, r);
			if (bin_prove_local(&pi, &ctx_pi, &com, m[SLOT_S], m[SLOT_F],
					wraw, g, &key, r, rng)) {
				break;
			}
		}
		/* Coordinate 0 of the randomness is the one the Ajtai row reaches and
		 * the message rows do not: b2 is zero there while B1 is one. Adding a
		 * single unit there leaves the rows and the norm bound untouched and
		 * breaks only the equation that ties z_mask to the commitment. */
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_set_coeff_ui(tst_zm[0][k], 0,
					nmod_add(nmod_poly_get_coeff_ui(tst_zm[0][k], 0), 1,
					tst_zm[0][k]->mod));
		}
		TEST_ASSERT(bin_verify_local(&pi, &ctx_pi, &com, &key) == 0, end);
	} TEST_END;

	TEST_ONCE("is_bin rejects a mask opening that is not short") {
		nmod_poly_zero(s);
		nmod_poly_fit_length(s, DEGREE);
		for (int i = 0; i < DEGREE; i++) {
			nmod_poly_set_coeff_ui(s, i, i & 1);
		}
		for (tries = 0; tries < 64; tries++) {
			for (int i = 0; i < LNP_WIDTH; i++) {
				commit_sample_short_crt(r[i]);
			}
			for (int i = 0; i < SLOTS; i++) {
				for (int k = 0; k < NCRT; k++) {
					nmod_poly_zero(m[i][k]);
				}
			}
			pcrt_poly_reduce(m[SLOT_S][0], s, 0);
			pcrt_poly_reduce(m[SLOT_S][1], s, 1);
			lnp_isbin_product(m[SLOT_F], m[SLOT_S]);
			lnp_sample_proj_mask(w, wraw);
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_set(m[SLOT_W][k], w[k]);
			}
			for (int i = 0; i < LNP_LAMBDA; i++) {
				lnp_sample_ct_zero(g[i], rng);
			}
			lnp_commit(&com, m, &key, r);
			if (bin_prove_local(&pi, &ctx_pi, &com, m[SLOT_S], m[SLOT_F],
					wraw, g, &key, r, rng)) {
				break;
			}
		}
		/* An unbounded opening is what solving the rows for z_mask produces,
		 * so the norm bound is what makes that solve useless. Both checks
		 * reject this input; isolating the norm bound the way the test above
		 * isolates the commitment equation would mean exhibiting a long
		 * vector in the kernel of the whole key, which is the MSIS problem
		 * the binding rests on. */
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_scalar_mul_nmod(tst_zm[0][k], tst_zm[0][k],
					(ulong) 1 << 20);
		}
		TEST_ASSERT(bin_verify_local(&pi, &ctx_pi, &com, &key) == 0, end);
	} TEST_END;

  end:
	nmod_poly_clear(s);
	nmod_poly_clear(wraw);
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
	for (int i = 0; i < LNP_LAMBDA; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(g[i][k]);
		}
	}
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(w[k]);
	}
	lnp_keyfree(&key);
	lnp_com_free(&com);
	lnp_binproof_free(&pi);
	lnp_binctx_free(&ctx_pi);
}

int main(int argc, char *argv[]) {
	flint_rand_t rand;

	flint_rand_init(rand);
	commit_setup();

	if (phase_selected(argc, argv, "test")) {
		printf("\n** Tests for the LNP proof machinery:\n\n");
		test_auto(rand);
			test_isbin(rand);
		test_adaptive(rand);
		test_isbin_full(rand);
		test_range(rand);
	}

	commit_finish();
	flint_rand_clear(rand);
	return 0;
}

#endif /* LNP_MAIN */
