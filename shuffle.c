#include <math.h>
#include <stdlib.h>

#include "param.h"
#include "commit.h"
#include "lnp.h"
#include "test.h"
#include "bench.h"
#include "assert.h"
#include "fastrandombytes.h"
#include "sha.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/* Number of messages to shuffle. Override at build time with -DMSGS=<n>; it
 * sizes stack arrays, so it cannot be a runtime parameter as things stand. */
#ifndef MSGS
#define MSGS        25
#endif

/*
 * The proof of shuffle below follows Neff's paradigm, but it does *not* use the
 * plain product identity
 *
 *     \prod (a_i - X) = \prod (b_i - X)                                    (*)
 *
 * that the original version of this code used. Over R_p = Z_p[x]/(x^N + 1) the
 * polynomial x^N + 1 always factors, so R_p is never a field and (*) does not
 * imply that the two lists are related by a permutation: it only implies that
 * they are permuted inside each CRT component, possibly by *different*
 * permutations. Bootle, Lyubashevsky and Merino-Gallardo, "Efficient Verifiable
 * Mixnets from Lattices, Revisited" (ePrint 2025/658), exploit exactly that gap
 * to break soundness.
 *
 * We therefore use the product from their Lemma 5, which ties every message to
 * its index and forces the per-component permutations to agree:
 *
 *     \prod (a_i + g(i) * X1 - X2) = \prod (b_i + sigma_i * X1 - X2),
 *
 * where g : [N] -> D is injective, D is a set whose pairwise differences are
 * invertible, and sigma_i \in D is committed by the prover before the
 * challenges X1, X2 are drawn. An honest prover sets sigma_i = g(pi(i)) for the
 * secret permutation pi.
 *
 * Lemma 5 additionally requires the committed sigma_i to lie in D. Protocol 1
 * of that paper discharges this with a sub-proof of is_bin(sigma_i), delegated
 * to a general-purpose proof system (LaZer). We instead prove that sigma_i is
 * *short*, which is what the lemma actually needs: for D we take a ball of
 * small norm, whose pairwise differences are invertible by Lemma 1, and we let
 * g map indices to monomials. The linear proof below therefore also masks the
 * message committed in p and the verifier bounds the norm of the response.
 *
 * CAVEAT: this is a deviation from Protocol 1, not merely an implementation of
 * it. The proof is relaxed, so what the extractor obtains is not "sigma_i is in
 * D" but "(d - d') * sigma_i is short" for a challenge difference d - d'. The
 * step of Lemma 5 that lifts sigma_i = g(j) from one CRT component to the whole
 * ring has to be redone accordingly: if sigma_i = g(j) mod p_l, then
 * (d - d') * (sigma_i - g(j)) is short and vanishes mod p_l, so by Lemma 1 it
 * is zero, and since d - d' is invertible sigma_i = g(j) over R. That argument
 * is sound but it is a restatement of the lemma, which has not been reviewed.
 * See the note in README.md.
 */

/**
 * The map g : [N] -> D of Lemma 5, instantiated as the monomial map i -> x^i.
 *
 * D is the set of monomials. The difference x^a - x^b of two distinct elements
 * has l2-norm sqrt(2), far below MODP^(1/2) ~ 62501, so by Lemma 1 it is
 * invertible in R_p, as Lemma 5 requires. The map is injective for i < DEGREE,
 * so it supports up to DEGREE messages.
 *
 * Monomials are used in place of the ring elements with binary coefficients of
 * the paper because multiplication by a monomial is a signed rotation, so
 * ||d * x^i|| = ||d|| exactly. That keeps the Gaussian mask SIGMA_S narrow
 * enough for the norm extracted from the proof below to stay under MODP^(1/2),
 * which a binary encoding would only barely achieve, and only for small N.
 *
 * @param[out] out			- the resulting ring element.
 * @param[in] i				- the index to encode, must be below DEGREE.
 */
static void int_to_monomial(nmod_poly_t out, int i) {
	nmod_poly_zero(out);
	nmod_poly_fit_length(out, DEGREE);
	nmod_poly_set_coeff_ui(out, i, 1);
}

/**
 * Test if a ring element is invertible in R_p, by checking that it is non-zero
 * in both CRT components (Fact 1).
 *
 * @param[in] a				- the ring element to test.
 * @return 1 if the element is invertible, 0 otherwise.
 */
static int is_invertible(nmod_poly_t a) {
	nmod_poly_t t;
	int result;

	nmod_poly_init(t, MODP);
	pcrt_poly_reduce(t, a, 0);
	result = !nmod_poly_is_zero(t);
	pcrt_poly_reduce(t, a, 1);
	result &= !nmod_poly_is_zero(t);
	nmod_poly_clear(t);
	return result;
}

/**
 * Absorb a polynomial into a hash computation in a canonical way.
 *
 * The coefficients are absorbed one by one up to DEGREE, so that the digest
 * depends only on the value of the polynomial. Absorbing the raw coefficient
 * array instead makes the digest depend on the capacity FLINT happens to have
 * allocated, which is a function of how the polynomial was built and therefore
 * need not agree between the prover and the verifier.
 *
 * @param[in,out] sha		- the hash context.
 * @param[in] p				- the polynomial to absorb.
 */
static void hash_poly(SHA256Context *sha, nmod_poly_t p) {
	uint64_t buf[DEGREE];

	for (int i = 0; i < DEGREE; i++) {
		buf[i] = nmod_poly_get_coeff_ui(p, i);
	}
	SHA256Input(sha, (const uint8_t *)buf, sizeof(buf));
}

/* Everything the is_bin sub-proof contributes for one message: the product
 * relation, the constant coefficient being zero, and the norm bound. The three
 * together say that the committed permutation element is binary, which is the
 * set D that Lemma 5 needs. */
typedef struct _isbin_t {
	lnpbinproof_t all;
	lnpbinctx_t ctx;
} isbin_t;

/*
 * The linear proof now relates *three* commitments instead of two. It proves
 * knowledge of openings of x, p and _x such that
 *
 *     alpha * <msg of x> + gamma * <msg of p> + beta = <msg of _x>,
 *
 * for public alpha, gamma, beta. The extra commitment p is what carries the
 * committed permutation element sigma_i of Lemma 5: in the shuffle, the factor
 * b_i = _m_i + sigma_i * tau - rho is no longer public, and splits into the
 * public part (_m_i - rho), folded into beta, and the committed part sigma_i
 * scaled by the public tau, folded into gamma.
 */

/*
 * The prover in two phases, so that one challenge can cover every message.
 *
 * Sharing a challenge is what lets the is_bin argument share the linear
 * proof's opening, and once the challenge is shared across messages the
 * rejection test has to span them too: testing each separately would multiply
 * the abort probabilities, 2^-MSGS rather than about a half. So the masks are
 * drawn at SIGMA_B, wide enough for a term sqrt(MSGS) times longer than one
 * message's.
 */
/* One challenge for the whole batch. It absorbs every commitment and every
 * first message, so no message's transcript can be replayed against another's,
 * and so the single opening each message sends answers a challenge that
 * depends on all of them. */
static void batch_hash(nmod_poly_t d[2], commitkey_t *key, lnpkey_t *lkey,
		commit_t com[MSGS], lnpcom_t p[MSGS], commit_t dc[MSGS],
		nmod_poly_t t[MSGS][2], nmod_poly_t tp[MSGS][2],
		nmod_poly_t _t[MSGS][2], nmod_poly_t u[MSGS][2], isbin_t ib[MSGS],
		nmod_poly_t beta, lnpmaskcom_t *mcom, lnpbatch_t *batch) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	uint32_t buf;

	SHA256Reset(&sha);
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				hash_poly(&sha, key->B1[i][j][k]);
				if (i == 0) {
					hash_poly(&sha, key->b2[j][k]);
				}
			}
		}
		for (int j = 0; j < LNP_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				hash_poly(&sha, lkey->B1[i][j][k]);
				for (int q = 0; q < SLOTS; q++) {
					hash_poly(&sha, lkey->b2[q][j][k]);
				}
			}
		}
	}
	hash_poly(&sha, beta);
	for (int k = 0; k < NCRT; k++) {
		for (int i = 0; i < HEIGHT; i++) {
			hash_poly(&sha, mcom->c1[i][k]);
			hash_poly(&sha, batch->w[i][k]);
		}
		for (int i = 0; i < LNP_LAMBDA; i++) {
			hash_poly(&sha, mcom->c2[i][k]);
			hash_poly(&sha, batch->h[i][k]);
			hash_poly(&sha, batch->v[i][k]);
		}
	}
	for (int l = 0; l < MSGS; l++) {
		for (int k = 0; k < NCRT; k++) {
			hash_poly(&sha, com[l].c1[k]);
			hash_poly(&sha, com[l].c2[k]);
			hash_poly(&sha, dc[l].c1[k]);
			hash_poly(&sha, dc[l].c2[k]);
			for (int i = 0; i < HEIGHT; i++) {
				hash_poly(&sha, p[l].c1[i][k]);
			}
			for (int i = 0; i < SLOTS; i++) {
				hash_poly(&sha, p[l].c2[i][k]);
			}
			hash_poly(&sha, t[l][k]);
			hash_poly(&sha, tp[l][k]);
			hash_poly(&sha, _t[l][k]);
			hash_poly(&sha, u[l][k]);
			hash_poly(&sha, ib[l].all.t[k]);
		}
		SHA256Input(&sha, (const uint8_t *)ib[l].all.zp,
				PROJ * sizeof(ulong));
	}
	SHA256Result(&sha, hash);

	fastrandombytes_setseed(hash);
	for (int i = 0; i < 2; i++) {
		nmod_poly_zero(d[i]);
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
	pcrt_poly_reduce(d[0], d[1], 0);
	pcrt_poly_reduce(d[1], d[1], 1);
}

static void lin_first(nmod_poly_t y[WIDTH][2], nmod_poly_t w[LNP_WIDTH][2],
		nmod_poly_t _y[WIDTH][2], nmod_poly_t t[2], nmod_poly_t tp[2],
		nmod_poly_t _t[2], nmod_poly_t u[2], lnpcom_t *p,
		commitkey_t *key, lnpkey_t *lkey, nmod_poly_t alpha,
		nmod_poly_t gamma, nmod_poly_t s[LNP_WIDTH][2], isbin_t *ib,
		nmod_poly_t sig[2], lnpbatch_t *batch) {
	nmod_poly_t tmp, a[2], g[2];

	nmod_poly_init(tmp, MODP);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(a[i], MODP);
		nmod_poly_init(g[i], MODP);
		pcrt_poly_reduce(a[i], alpha, i);
		pcrt_poly_reduce(g[i], gamma, i);
		nmod_poly_zero(t[i]);
		nmod_poly_zero(tp[i]);
		nmod_poly_zero(_t[i]);
		nmod_poly_zero(u[i]);
	}

	for (int i = 0; i < WIDTH; i++) {
		commit_sample_gauss_batch_crt(y[i]);
		commit_sample_gauss_batch_crt(_y[i]);
	}
	for (int i = 0; i < LNP_WIDTH; i++) {
		commit_sample_gauss_batch_crt(w[i]);
	}
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_mulmod(tmp, key->B1[i][j][k], y[j][k], k);
				nmod_poly_add(t[k], t[k], tmp);
				pcrt_poly_mulmod(tmp, key->B1[i][j][k], _y[j][k], k);
				nmod_poly_add(_t[k], _t[k], tmp);
			}
		}
		for (int j = 0; j < LNP_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_mulmod(tmp, lkey->B1[i][j][k], w[j][k], k);
				nmod_poly_add(tp[k], tp[k], tmp);
			}
		}
	}

	/* u = alpha * <b2, y> + gamma * <b2_lnp, w> - <b2, y'>. */
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			pcrt_poly_mulmod(tmp, key->b2[i][j], y[i][j], j);
			pcrt_poly_mulmod(tmp, tmp, a[j], j);
			nmod_poly_add(u[j], u[j], tmp);
			pcrt_poly_mulmod(tmp, key->b2[i][j], _y[i][j], j);
			nmod_poly_sub(u[j], u[j], tmp);
		}
	}
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			pcrt_poly_mulmod(tmp, lkey->b2[SLOT_S][i][j], w[i][j], j);
			pcrt_poly_mulmod(tmp, tmp, g[j], j);
			nmod_poly_add(u[j], u[j], tmp);
		}
	}

	/* The is_bin first messages ride on the same mask. */
	lnp_bin_first(&ib->all, &ib->ctx, batch, p, sig, lkey, s, w);

	nmod_poly_clear(tmp);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(a[i]);
		nmod_poly_clear(g[i]);
	}
}

/* Answer the shared challenge, accumulating what the one rejection test over
 * every message's responses will need. */
static void lin_respond(nmod_poly_t y[WIDTH][2], nmod_poly_t w[LNP_WIDTH][2],
		nmod_poly_t _y[WIDTH][2], nmod_poly_t d[2], nmod_poly_t r[WIDTH][2],
		nmod_poly_t s[LNP_WIDTH][2], nmod_poly_t _r[WIDTH][2],
		int64_t *dot, int64_t *norm) {
	nmod_poly_t dr[WIDTH][2], ds[LNP_WIDTH][2], _dr[WIDTH][2];

	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_init(dr[i][j], MODP);
			nmod_poly_init(_dr[i][j], MODP);
			pcrt_poly_mulmod(dr[i][j], d[j], r[i][j], j);
			nmod_poly_add(y[i][j], y[i][j], dr[i][j]);
			pcrt_poly_mulmod(_dr[i][j], d[j], _r[i][j], j);
			nmod_poly_add(_y[i][j], _y[i][j], _dr[i][j]);
		}
	}
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_init(ds[i][j], MODP);
			pcrt_poly_mulmod(ds[i][j], d[j], s[i][j], j);
			nmod_poly_add(w[i][j], w[i][j], ds[i][j]);
		}
	}
	commit_rej_accumulate(dot, norm, y, dr, WIDTH);
	commit_rej_accumulate(dot, norm, _y, _dr, WIDTH);
	commit_rej_accumulate(dot, norm, w, ds, LNP_WIDTH);

	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_clear(dr[i][j]);
			nmod_poly_clear(_dr[i][j]);
		}
	}
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_clear(ds[i][j]);
		}
	}
}

static int lin_verifier(nmod_poly_t y[WIDTH][2], nmod_poly_t w[LNP_WIDTH][2],
		nmod_poly_t _y[WIDTH][2], nmod_poly_t t[2],
		nmod_poly_t tp[2], nmod_poly_t _t[2],
		nmod_poly_t u[2], commit_t x, lnpcom_t *p, commit_t _x,
		commitkey_t *key, lnpkey_t *lkey, nmod_poly_t alpha,
		nmod_poly_t gamma, nmod_poly_t beta, isbin_t *ib,
		nmod_poly_t _d[2], lnpmaskcom_t *mcom, pcrt_poly_t acc[LNP_LAMBDA]) {
	nmod_poly_t tmp, a[2], g[2], b[2];
	nmod_poly_t v[2], vp[2], _v[2], lhs[2], rhs[2];
	nmod_poly_t z[WIDTH], zp[LNP_WIDTH], _z[WIDTH];
	int result = 1;

	nmod_poly_init(tmp, MODP);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(z[i], MODP);
		nmod_poly_init(_z[i], MODP);
	}
	for (int i = 0; i < LNP_WIDTH; i++) {
		nmod_poly_init(zp[i], MODP);
	}
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(a[i], MODP);
		nmod_poly_init(g[i], MODP);
		nmod_poly_init(b[i], MODP);
		nmod_poly_init(v[i], MODP);
		nmod_poly_init(vp[i], MODP);
		nmod_poly_init(_v[i], MODP);
		nmod_poly_init(lhs[i], MODP);
		nmod_poly_init(rhs[i], MODP);
		nmod_poly_zero(v[i]);
		nmod_poly_zero(vp[i]);
		nmod_poly_zero(_v[i]);
		/* Reduce the public coefficients into CRT representation once. */
		pcrt_poly_reduce(a[i], alpha, i);
		pcrt_poly_reduce(g[i], gamma, i);
		pcrt_poly_reduce(b[i], beta, i);
	}

	/* Sample challenge. */

	/* The same opening w answers the is_bin argument, whose aggregated half is
	 * settled once for the batch. */
	result &= lnp_bin_check(&ib->all, &ib->ctx, p, lkey, _d, w, acc);

	/* Verifier checks norms, reconstructing from CRT representation. These are
	 * soundness checks and must make verification fail, so they are ordinary
	 * checks rather than assertions that vanish under NDEBUG. */
	for (int i = 0; i < WIDTH; i++) {
		pcrt_poly_rec(z[i], y[i]);
		pcrt_poly_rec(_z[i], _y[i]);
		result &= commit_norm2_leq(z[i], (uint64_t) 4 * DEGREE * SIGMA_B * SIGMA_B);
		result &= commit_norm2_leq(_z[i], (uint64_t) 4 * DEGREE * SIGMA_B * SIGMA_B);
	}
	for (int i = 0; i < LNP_WIDTH; i++) {
		pcrt_poly_rec(zp[i], w[i]);
		result &= commit_norm2_leq(zp[i], (uint64_t) 4 * DEGREE * SIGMA_B * SIGMA_B);
	}

	/* The set D that Lemma 5 requires is now the binary ring elements, and
	 * membership is established by the is_bin proof that shuffle_verifier runs
	 * against this same commitment, not by a norm check here. */

	/* Verifier computes B1z, B1z_p and B1z'. */
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_mulmod(tmp, key->B1[i][j][k], y[j][k], k);
				nmod_poly_add(v[k], v[k], tmp);
				pcrt_poly_mulmod(tmp, key->B1[i][j][k], _y[j][k], k);
				nmod_poly_add(_v[k], _v[k], tmp);
			}
		}
		for (int j = 0; j < LNP_WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_mulmod(tmp, lkey->B1[i][j][k], w[j][k], k);
				nmod_poly_add(vp[k], vp[k], tmp);
			}
		}
	}
	/* Verifier checks that B1z = t + d[x], B1z_p = t_p + d[p], B1z' = t' + d[x']. */
	for (int j = 0; j < NCRT; j++) {
		pcrt_poly_mulmod(tmp, _d[j], x.c1[j], j);
		nmod_poly_add(lhs[j], t[j], tmp);
		result &= nmod_poly_equal(lhs[j], v[j]);
		pcrt_poly_mulmod(tmp, _d[j], p->c1[0][j], j);
		nmod_poly_add(lhs[j], tp[j], tmp);
		result &= nmod_poly_equal(lhs[j], vp[j]);
		pcrt_poly_mulmod(tmp, _d[j], _x.c1[j], j);
		nmod_poly_add(lhs[j], _t[j], tmp);
		result &= nmod_poly_equal(lhs[j], _v[j]);
	}

	/* Verifier checks the linear relation
	 *     alpha * <b2, z> + gamma * <b2, z_p> - <b2, z'>
	 *         = u + d * (alpha * x.c2 + gamma * p.c2 - x'.c2 + beta). */
	for (int j = 0; j < NCRT; j++) {
		pcrt_poly_mulmod(lhs[j], a[j], x.c2[j], j);
		pcrt_poly_mulmod(tmp, g[j], p->c2[SLOT_S][j], j);
		nmod_poly_add(lhs[j], lhs[j], tmp);
		nmod_poly_sub(lhs[j], lhs[j], _x.c2[j]);
		nmod_poly_add(lhs[j], lhs[j], b[j]);
		pcrt_poly_mulmod(lhs[j], lhs[j], _d[j], j);
		nmod_poly_add(lhs[j], lhs[j], u[j]);
		nmod_poly_zero(rhs[j]);
	}

	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			pcrt_poly_mulmod(tmp, key->b2[i][j], y[i][j], j);
			pcrt_poly_mulmod(tmp, a[j], tmp, j);
			nmod_poly_add(rhs[j], rhs[j], tmp);
			pcrt_poly_mulmod(tmp, key->b2[i][j], _y[i][j], j);
			nmod_poly_sub(rhs[j], rhs[j], tmp);
		}
	}
	for (int i = 0; i < LNP_WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			pcrt_poly_mulmod(tmp, lkey->b2[SLOT_S][i][j], w[i][j], j);
			pcrt_poly_mulmod(tmp, g[j], tmp, j);
			nmod_poly_add(rhs[j], rhs[j], tmp);
		}
	}
	for (int j = 0; j < NCRT; j++) {
		result &= nmod_poly_equal(lhs[j], rhs[j]);
	}

	nmod_poly_clear(tmp);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_clear(z[i]);
		nmod_poly_clear(_z[i]);
	}
	for (int i = 0; i < LNP_WIDTH; i++) {
		nmod_poly_clear(zp[i]);
	}
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(a[i]);
		nmod_poly_clear(g[i]);
		nmod_poly_clear(b[i]);
		nmod_poly_clear(v[i]);
		nmod_poly_clear(vp[i]);
		nmod_poly_clear(_v[i]);
		nmod_poly_clear(lhs[i]);
		nmod_poly_clear(rhs[i]);
	}
	return result;
}

void shuffle_hash(nmod_poly_t beta, commit_t c[MSGS], lnpcom_t p[MSGS],
		commit_t d[MSGS], nmod_poly_t _m[MSGS], nmod_poly_t tau,
		nmod_poly_t rho) {
	flint_rand_t rand;
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	uint64_t seed0, seed1, seed2, seed3;

	SHA256Reset(&sha);

	for (int i = 0; i < MSGS; i++) {
		hash_poly(&sha, _m[i]);
		for (int j = 0; j < NCRT; j++) {
			hash_poly(&sha, c[i].c1[j]);
			hash_poly(&sha, c[i].c2[j]);
			hash_poly(&sha, p[i].c1[0][j]);
			hash_poly(&sha, p[i].c2[SLOT_S][j]);
			hash_poly(&sha, d[i].c1[j]);
			hash_poly(&sha, d[i].c2[j]);
		}
	}
	hash_poly(&sha, tau);
	hash_poly(&sha, rho);
	SHA256Result(&sha, hash);

	flint_rand_init(rand);
	memcpy(&seed0, hash, sizeof(uint64_t));
	memcpy(&seed1, hash + sizeof(uint64_t), sizeof(uint64_t));
	memcpy(&seed2, hash + 2 * sizeof(uint64_t), sizeof(uint64_t));
	memcpy(&seed3, hash + 3 * sizeof(uint64_t), sizeof(uint64_t));
	seed0 ^= seed2;
	seed1 ^= seed3;
	flint_rand_set_seed(rand, seed0, seed1);
	/* Challenges are sampled of degree below DEGCRT = N/2, so that the
	 * difference of two distinct challenges is non-zero in both CRT components
	 * and hence invertible. */
	commit_sample_rand(beta, rand, DEGCRT);
	flint_rand_clear(rand);
}

/**
 * Compute the public coefficients of the l-th linear relation of the product
 * argument, namely
 *
 *     alpha * a_l + gamma_raw * b_l = <message committed in D_l>,
 *
 * with a_l = m_l + g(l) * tau - rho and b_l = _m_l + sigma_l * tau - rho.
 * Only m_l (inside com[l]) and sigma_l (inside P[l]) are secret, so the
 * relation handed to the linear proof is
 *
 *     alpha * m_l + (gamma_raw * tau) * sigma_l + pub = <message of D_l>,
 *     pub = alpha * (g(l) * tau - rho) + gamma_raw * (_m_l - rho).
 *
 * The (-1)^N sign of the last equation is folded into gamma_raw, so that prover
 * and verifier treat every index uniformly.
 *
 * @param[out] alpha		- the coefficient of the message of com[l].
 * @param[out] gamma		- the coefficient of the message of P[l].
 * @param[out] pub			- the public additive term.
 * @param[in] l				- the index of the relation.
 * @param[in] s				- the values s_i sent by the prover.
 * @param[in] _m			- the public output list of messages.
 * @param[in] beta			- the challenge of the product argument.
 * @param[in] tau			- the challenge X1 of Lemma 5.
 * @param[in] rho			- the challenge X2 of Lemma 5.
 */
static void shuffle_coeffs(nmod_poly_t alpha, nmod_poly_t gamma,
		nmod_poly_t pub, int l, nmod_poly_t s[MSGS], nmod_poly_t _m[MSGS],
		nmod_poly_t beta, nmod_poly_t tau, nmod_poly_t rho) {
	nmod_poly_t raw, t0, t1;

	nmod_poly_init(raw, MODP);
	nmod_poly_init(t0, MODP);
	nmod_poly_init(t1, MODP);

	if (l == 0) {
		nmod_poly_set(alpha, beta);
	} else {
		nmod_poly_set(alpha, s[l - 1]);
	}

	if (l < MSGS - 1) {
		nmod_poly_set(raw, s[l]);
	} else {
		/* Coefficient (-1)^N of b_{N-1} in the last equation. */
		if (MSGS & 1) {
			nmod_poly_zero(raw);
			nmod_poly_sub(raw, raw, beta);
		} else {
			nmod_poly_set(raw, beta);
		}
	}

	/* pub = alpha * (g(l) * tau - rho) + raw * (_m_l - rho). */
	int_to_monomial(t0, l);
	commit_poly_mulmod(t0, t0, tau);
	nmod_poly_sub(t0, t0, rho);
	commit_poly_mulmod(pub, alpha, t0);
	nmod_poly_sub(t1, _m[l], rho);
	commit_poly_mulmod(t1, raw, t1);
	nmod_poly_add(pub, pub, t1);

	/* The committed sigma_l enters scaled by tau. */
	commit_poly_mulmod(gamma, raw, tau);

	nmod_poly_clear(raw);
	nmod_poly_clear(t0);
	nmod_poly_clear(t1);
}

static void shuffle_prover(nmod_poly_t y[MSGS][WIDTH][2],
		nmod_poly_t w[MSGS][LNP_WIDTH][2], nmod_poly_t _y[MSGS][WIDTH][2],
		nmod_poly_t t[MSGS][2],
		nmod_poly_t tp[MSGS][2], nmod_poly_t _t[MSGS][2],
		nmod_poly_t u[MSGS][2], commit_t d[MSGS],
		nmod_poly_t s[MSGS], commit_t com[MSGS], lnpcom_t p[MSGS],
		isbin_t ib[MSGS], nmod_poly_t m[MSGS], nmod_poly_t _m[MSGS],
		nmod_poly_t sigma[MSGS], nmod_poly_t r[MSGS][WIDTH][2],
		nmod_poly_t pr[MSGS][LNP_WIDTH][2], nmod_poly_t tau, nmod_poly_t rho,
		commitkey_t *key, lnpkey_t *lkey, lnpmaskkey_t *mkey,
		lnpmaskcom_t *mcom, lnpbatch_t *batch,
		nmod_poly_t mr[MASK_WIDTH][2], nmod_poly_t zm[MASK_WIDTH][2],
		flint_rand_t rng) {
	nmod_poly_t beta, alpha, gamma, pub, t0, t1;
	nmod_poly_t a[MSGS], b[MSGS], theta[MSGS], _r[MSGS][WIDTH][2];
	nmod_poly_t sig[MSGS][2];

	nmod_poly_init(t0, MODP);
	nmod_poly_init(t1, MODP);
	nmod_poly_init(beta, MODP);
	nmod_poly_init(alpha, MODP);
	nmod_poly_init(gamma, MODP);
	nmod_poly_init(pub, MODP);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(a[i], MODP);
		nmod_poly_init(b[i], MODP);
		nmod_poly_init(theta[i], MODP);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(sig[i][k], MODP);
			/* CRT representation of the committed permutation element. */
			pcrt_poly_reduce(sig[i][k], sigma[i], k);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_init(_r[i][j][k], MODP);
			}
		}
	}

	/* Build the factors of the product from Lemma 5:
	 *     a_i = m_i + g(i) * tau - rho,
	 *     b_i = _m_i + sigma_i * tau - rho,
	 * with sigma_i = g(pi(i)) committed in P[i]. Tying the index to every
	 * message is what forces the permutations of the two CRT components to
	 * coincide, and is what the original product a_i - rho lacked. */
	for (int i = 0; i < MSGS; i++) {
		int_to_monomial(t1, i);
		commit_poly_mulmod(t0, t1, tau);
		nmod_poly_add(a[i], m[i], t0);
		nmod_poly_sub(a[i], a[i], rho);

		commit_poly_mulmod(t0, sigma[i], tau);
		nmod_poly_add(b[i], _m[i], t0);
		nmod_poly_sub(b[i], b[i], rho);
		/* The product argument inverts the b_i. This fails only with the
		 * negligible probability accounted for in the completeness error. */
		assert(is_invertible(b[i]));
	}

	/* Prover samples theta_i and computes commitments D_i. */
	commit_sample_rand(theta[0], rng, DEGREE);
	commit_poly_mulmod(t0, theta[0], b[0]);
	for (int j = 0; j < WIDTH; j++) {
		commit_sample_short_crt(_r[0][j]);
	}
	commit_doit(&d[0], t0, key, _r[0]);
	for (int i = 1; i < MSGS - 1; i++) {
		commit_sample_rand(theta[i], rng, DEGREE);
		commit_poly_mulmod(t0, theta[i - 1], a[i]);
		commit_poly_mulmod(t1, theta[i], b[i]);
		nmod_poly_add(t0, t0, t1);
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short_crt(_r[i][j]);
		}
		commit_doit(&d[i], t0, key, _r[i]);
	}
	commit_poly_mulmod(t0, theta[MSGS - 2], a[MSGS - 1]);
	for (int j = 0; j < WIDTH; j++) {
		commit_sample_short_crt(_r[MSGS - 1][j]);
	}
	commit_doit(&d[MSGS - 1], t0, key, _r[MSGS - 1]);

	shuffle_hash(beta, com, p, d, _m, tau, rho);
	commit_poly_mulmod(s[0], theta[0], b[0]);
	commit_poly_mulmod(t0, beta, a[0]);
	nmod_poly_sub(s[0], s[0], t0);
	nmod_poly_invmod(t0, b[0], *commit_poly());
	commit_poly_mulmod(s[0], s[0], t0);
	for (int i = 1; i < MSGS - 1; i++) {
		commit_poly_mulmod(s[i], theta[i - 1], a[i]);
		commit_poly_mulmod(t0, theta[i], b[i]);
		nmod_poly_add(s[i], s[i], t0);
		commit_poly_mulmod(t0, s[i - 1], a[i]);
		nmod_poly_sub(s[i], s[i], t0);
		nmod_poly_invmod(t0, b[i], *commit_poly());
		commit_poly_mulmod(s[i], s[i], t0);
	}

	/* The setup half of the is_bin argument, once per message. It has to run
	 * before the linear proof, which drives the mask-dependent half and needs
	 * the context this produces. */
	{
		pcrt_poly_t msg[SLOTS], g[LNP_LAMBDA], wc;
		nmod_poly_t wraw;

		nmod_poly_init(wraw, MODP);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(wc[k], MODP);
		}
		for (int i = 0; i < SLOTS; i++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(msg[i][k], MODP);
			}
		}
		for (int i = 0; i < LNP_LAMBDA; i++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(g[i][k], MODP);
			}
		}
		/* One set of masks covers the whole batch, so they are committed once
		 * and their contribution to the aggregated values added once, before
		 * any message adds its own share. */
		for (int i = 0; i < LNP_LAMBDA; i++) {
			lnp_sample_ct_zero(g[i], rng);
		}
		lnp_mask_commit(mcom, g, mkey, mr);
		lnp_batch_zero(batch);
		for (int j = 0; j < LNP_LAMBDA; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_set(batch->h[j][k], g[j][k]);
			}
		}

		for (int l = 0; l < MSGS; l++) {
			/* The setup half of the is_bin argument fixes the projection,
			 * which cannot be retried once its mask is committed, so a
			 * rejection here means recommitting. The mask-dependent half runs
			 * inside lin_first, which owns the shared opening. */
			for (int tries = 0; tries < 64; tries++) {
				for (int i = 0; i < SLOTS; i++) {
					for (int k = 0; k < NCRT; k++) {
						nmod_poly_zero(msg[i][k]);
					}
				}
				for (int k = 0; k < NCRT; k++) {
					nmod_poly_set(msg[SLOT_S][k], sig[l][k]);
				}
				lnp_isbin_product(msg[SLOT_F], msg[SLOT_S]);
				lnp_sample_proj_mask(wc, wraw);
				for (int k = 0; k < NCRT; k++) {
					nmod_poly_set(msg[SLOT_W][k], wc[k]);
				}
				lnp_commit(&p[l], msg, lkey, pr[l]);
				if (lnp_bin_setup(&ib[l].all, &ib[l].ctx, batch, &p[l], mcom,
						msg[SLOT_S], msg[SLOT_F], wraw, lkey)) {
					break;
				}
			}
		}
		nmod_poly_clear(wraw);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(wc[k]);
		}
		for (int i = 0; i < SLOTS; i++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(msg[i][k]);
			}
		}
		for (int i = 0; i < LNP_LAMBDA; i++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(g[i][k]);
			}
		}
	}

	/* One challenge for every message, so one rejection test over all their
	 * responses together. */
	{
		nmod_poly_t dch[NCRT];
		int64_t dot, norm;
		uint64_t sigma_sqr = (uint64_t) SIGMA_B * SIGMA_B;
		int rej;

		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(dch[k], MODP);
		}
		nmod_poly_t ym[MASK_WIDTH][2], cm[MASK_WIDTH][2];

		for (int i = 0; i < MASK_WIDTH; i++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(ym[i][k], MODP);
				nmod_poly_init(cm[i][k], MODP);
			}
		}
		do {
			dot = norm = 0;
			for (int j = 0; j < LNP_LAMBDA; j++) {
				for (int k = 0; k < NCRT; k++) {
					nmod_poly_zero(batch->v[j][k]);
				}
			}
			for (int i = 0; i < MASK_WIDTH; i++) {
				commit_sample_gauss_batch_crt(ym[i]);
			}
			for (int l = 0; l < MSGS; l++) {
				shuffle_coeffs(alpha, gamma, pub, l, s, _m, beta, tau, rho);
				lin_first(y[l], w[l], _y[l], t[l], tp[l], _t[l], u[l],
						&p[l], key, lkey, alpha, gamma, pr[l], &ib[l],
						sig[l], batch);
			}
			lnp_batch_first(batch, mkey, ym);
			batch_hash(dch, key, lkey, com, p, d, t, tp, _t, u, ib, beta,
					mcom, batch);
			for (int l = 0; l < MSGS; l++) {
				lin_respond(y[l], w[l], _y[l], dch, r[l], pr[l], _r[l],
						&dot, &norm);
			}
			for (int i = 0; i < MASK_WIDTH; i++) {
				for (int k = 0; k < NCRT; k++) {
					pcrt_poly_mulmod(cm[i][k], dch[k], mr[i][k], k);
					nmod_poly_add(zm[i][k], ym[i][k], cm[i][k]);
				}
			}
			commit_rej_accumulate(&dot, &norm, zm, cm, MASK_WIDTH);
			rej = commit_rej_decide(dot, norm, sigma_sqr);
		} while (rej);
		for (int i = 0; i < MASK_WIDTH; i++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(ym[i][k]);
				nmod_poly_clear(cm[i][k]);
			}
		}
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(dch[k]);
		}
	}


	nmod_poly_clear(t0);
	nmod_poly_clear(t1);
	nmod_poly_clear(beta);
	nmod_poly_clear(alpha);
	nmod_poly_clear(gamma);
	nmod_poly_clear(pub);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_clear(a[i]);
		nmod_poly_clear(b[i]);
		nmod_poly_clear(theta[i]);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(sig[i][k]);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_clear(_r[i][j][k]);
			}
		}
	}
}

static int shuffle_verifier(nmod_poly_t y[MSGS][WIDTH][2],
		nmod_poly_t w[MSGS][LNP_WIDTH][2], nmod_poly_t _y[MSGS][WIDTH][2],
		nmod_poly_t t[MSGS][2],
		nmod_poly_t tp[MSGS][2], nmod_poly_t _t[MSGS][2],
		nmod_poly_t u[MSGS][2], commit_t d[MSGS],
		nmod_poly_t s[MSGS], commit_t com[MSGS], lnpcom_t p[MSGS],
		isbin_t ib[MSGS], nmod_poly_t _m[MSGS], nmod_poly_t tau,
		nmod_poly_t rho, commitkey_t *key, lnpkey_t *lkey,
		lnpmaskkey_t *mkey, lnpmaskcom_t *mcom, lnpbatch_t *batch,
		nmod_poly_t zm[MASK_WIDTH][2]) {
	int result = 1;
	nmod_poly_t beta, alpha, gamma, pub, dch[NCRT];
	pcrt_poly_t acc[LNP_LAMBDA];

	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(dch[k], MODP);
	}
	for (int j = 0; j < LNP_LAMBDA; j++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(acc[j][k], MODP);
		}
	}
	nmod_poly_init(beta, MODP);
	nmod_poly_init(alpha, MODP);
	nmod_poly_init(gamma, MODP);
	nmod_poly_init(pub, MODP);

	shuffle_hash(beta, com, p, d, _m, tau, rho);
	for (int l = 0; l < MSGS; l++) {
		lnp_bin_public(&ib[l].ctx, &ib[l].all, &p[l], mcom, lkey);
	}
	batch_hash(dch, key, lkey, com, p, d, t, tp, _t, u, ib, beta, mcom,
			batch);
	/* Now verify each \Prod_LIN instance, one for each commitment. */
	for (int l = 0; l < MSGS; l++) {
		shuffle_coeffs(alpha, gamma, pub, l, s, _m, beta, tau, rho);
		/* This also checks membership of sigma_l in the set D that Lemma 5
		 * requires, since the is_bin argument now shares the linear proof's
		 * challenge and opening. */
		result &=
				lin_verifier(y[l], w[l], _y[l], t[l], tp[l], _t[l],
				u[l], com[l], &p[l], d[l], key, lkey, alpha, gamma, pub,
				&ib[l], dch, mcom, acc);
	}

	result &= lnp_batch_check(batch, mcom, mkey, dch, zm, acc);

	nmod_poly_clear(beta);
	nmod_poly_clear(alpha);
	nmod_poly_clear(gamma);
	nmod_poly_clear(pub);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(dch[k]);
	}
	for (int j = 0; j < LNP_LAMBDA; j++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(acc[j][k]);
		}
	}
	return result;
}

/**
 * Run a full proof of shuffle.
 *
 * @param[in] com			- the commitments to the input messages.
 * @param[in] m				- the input messages.
 * @param[in] _m			- the public output (shuffled) messages.
 * @param[in] sigma			- the permutation elements claimed by the prover. An
 *							  honest prover sets sigma[i] = g(pi(i)) for the
 *							  permutation with _m[i] = m[pi[i]].
 * @param[in] r				- the randomness of the input commitments.
 * @param[in] key			- the commitment key.
 * @param[in] rng			- the random number generator.
 * @return 1 if the proof verifies, 0 otherwise.
 */
static int run(commit_t com[MSGS], nmod_poly_t m[MSGS], nmod_poly_t _m[MSGS],
		nmod_poly_t sigma[MSGS], nmod_poly_t r[MSGS][WIDTH][2],
		commitkey_t *key, lnpkey_t *lkey, flint_rand_t rng) {
	int result = 1;
	commit_t d[MSGS];
	/* These are large enough that they do not belong on the stack. */
	static lnpcom_t p[MSGS];
	static isbin_t ib[MSGS];
	static nmod_poly_t pr[MSGS][LNP_WIDTH][2];
	static nmod_poly_t w[MSGS][LNP_WIDTH][2];
	static nmod_poly_t mr[MASK_WIDTH][2], zm[MASK_WIDTH][2];
	lnpmaskkey_t mkey;
	lnpmaskcom_t mcom;
	lnpbatch_t batch;
	nmod_poly_t tau, rho, s[MSGS], u[MSGS][2];
	nmod_poly_t y[MSGS][WIDTH][2], _y[MSGS][WIDTH][2];
	nmod_poly_t t[MSGS][2], tp[MSGS][2], _t[MSGS][2];
	pcrt_poly_t msg[SLOTS];

	nmod_poly_init(tau, MODP);
	nmod_poly_init(rho, MODP);
	lnp_maskkey_init(&mkey);
	lnp_maskkey_gen(&mkey, rng);
	lnp_maskcom_init(&mcom);
	lnp_batch_init(&batch);
	for (int i = 0; i < MASK_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(mr[i][k], MODP);
			nmod_poly_init(zm[i][k], MODP);
		}
		commit_sample_short_crt(mr[i]);
	}
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(msg[i][k], MODP);
			nmod_poly_zero(msg[i][k]);
		}
	}
	for (int i = 0; i < MSGS; i++) {
		commit_init(&d[i]);
		lnp_com_init(&p[i]);
		lnp_binproof_init(&ib[i].all);
		lnp_binctx_init(&ib[i].ctx);
		nmod_poly_init(s[i], MODP);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(t[i][k], MODP);
			nmod_poly_init(tp[i][k], MODP);
			nmod_poly_init(_t[i][k], MODP);
			nmod_poly_init(u[i][k], MODP);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_init(y[i][j][k], MODP);
				nmod_poly_init(_y[i][j][k], MODP);
			}
			for (int j = 0; j < LNP_WIDTH; j++) {
				nmod_poly_init(w[i][j][k], MODP);
				nmod_poly_init(pr[i][j][k], MODP);
			}
		}
	}

	/* Prover commits to the permutation elements sigma_i of Lemma 5. This must
	 * happen *before* the challenges tau and rho are drawn: the soundness
	 * argument needs the sigma_i to be fixed independently of them. The
	 * commitment is the LNP one, so that the is_bin sub-proof and the linear
	 * proof speak about the same object and no linking proof is needed. Only
	 * its Ajtai part and slot SLOT_S enter the shuffle hash, which leaves the
	 * proof material in the other slots free to be rewritten later. */
	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < LNP_WIDTH; j++) {
			commit_sample_short_crt(pr[i][j]);
		}
		pcrt_poly_reduce(msg[SLOT_S][0], sigma[i], 0);
		pcrt_poly_reduce(msg[SLOT_S][1], sigma[i], 1);
		lnp_commit(&p[i], msg, lkey, pr[i]);
	}

	/* Verifier samples the two challenges X1 = tau and X2 = rho of Lemma 5,
	 * of degree below N/2 so that differences of distinct challenges are
	 * invertible. */
	commit_sample_rand(tau, rng, DEGCRT);
	commit_sample_rand(rho, rng, DEGCRT);

	shuffle_prover(y, w, _y, t, tp, _t, u, d, s, com, p, ib, m, _m, sigma,
			r, pr, tau, rho, key, lkey, &mkey, &mcom, &batch, mr, zm, rng);

	result = shuffle_verifier(y, w, _y, t, tp, _t, u, d, s, com, p, ib, _m,
			tau, rho, key, lkey, &mkey, &mcom, &batch, zm);

	nmod_poly_clear(tau);
	nmod_poly_clear(rho);
	for (int i = 0; i < SLOTS; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(msg[i][k]);
		}
	}
	for (int i = 0; i < MASK_WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(mr[i][k]);
			nmod_poly_clear(zm[i][k]);
		}
	}
	lnp_maskkey_free(&mkey);
	lnp_maskcom_free(&mcom);
	lnp_batch_free(&batch);
	for (int i = 0; i < MSGS; i++) {
		commit_free(&d[i]);
		lnp_com_free(&p[i]);
		lnp_binproof_free(&ib[i].all);
		lnp_binctx_free(&ib[i].ctx);
		nmod_poly_clear(s[i]);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(t[i][k]);
			nmod_poly_clear(tp[i][k]);
			nmod_poly_clear(_t[i][k]);
			nmod_poly_clear(u[i][k]);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_clear(y[i][j][k]);
				nmod_poly_clear(_y[i][j][k]);
			}
			for (int j = 0; j < LNP_WIDTH; j++) {
				nmod_poly_clear(w[i][j][k]);
				nmod_poly_clear(pr[i][j][k]);
			}
		}
	}

	return result;
}

/**
 * Compute the product \prod (v_i - rho) that the original proof of shuffle
 * relied on, used by the tests to exhibit the attack.
 *
 * @param[out] out			- the resulting product.
 * @param[in] v				- the list of messages.
 * @param[in] rho			- the evaluation point.
 */
static void neff_product(nmod_poly_t out, nmod_poly_t v[MSGS],
		nmod_poly_t rho) {
	nmod_poly_t t, acc;

	nmod_poly_init(t, MODP);
	nmod_poly_init(acc, MODP);
	nmod_poly_zero(acc);
	nmod_poly_set_coeff_ui(acc, 0, 1);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_sub(t, v[i], rho);
		commit_poly_mulmod(acc, acc, t);
	}
	nmod_poly_set(out, acc);
	nmod_poly_clear(t);
	nmod_poly_clear(acc);
}

/**
 * Swap the first CRT component of two ring elements, leaving the second one
 * untouched. This is the manipulation of Section 4.1 of ePrint 2025/658: the
 * outputs are a permutation of the inputs inside each of the two fields the
 * ring splits into, but not over the ring itself.
 *
 * @param[out] out0			- the first manipulated element.
 * @param[out] out1			- the second manipulated element.
 * @param[in] in0			- the first input element.
 * @param[in] in1			- the second input element.
 */
static void crt_swap(nmod_poly_t out0, nmod_poly_t out1, nmod_poly_t in0,
		nmod_poly_t in1) {
	pcrt_poly_t a, b;

	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(a[i], MODP);
		nmod_poly_init(b[i], MODP);
	}
	/* out0 = (in1 mod p0, in0 mod p1), out1 = (in0 mod p0, in1 mod p1). */
	pcrt_poly_reduce(a[0], in1, 0);
	pcrt_poly_reduce(a[1], in0, 1);
	pcrt_poly_reduce(b[0], in0, 0);
	pcrt_poly_reduce(b[1], in1, 1);
	pcrt_poly_rec(out0, a);
	pcrt_poly_rec(out1, b);

	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(a[i]);
		nmod_poly_clear(b[i]);
	}
}

static void test(flint_rand_t rand) {
	commitkey_t key;
	lnpkey_t lkey;
	commit_t com[MSGS];
	int pi[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS], am[MSGS], sigma[MSGS], r[MSGS][WIDTH][2];
	nmod_poly_t rho, p0, p1, g0, g1;

	nmod_poly_init(rho, MODP);
	nmod_poly_init(p0, MODP);
	nmod_poly_init(p1, MODP);
	nmod_poly_init(g0, MODP);
	nmod_poly_init(g1, MODP);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
		nmod_poly_init(am[i], MODP);
		nmod_poly_init(sigma[i], MODP);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(r[i][j][k], MODP);
			}
		}
	}

	/* Generate commitment key-> */
	commit_setup();
	commit_keyinit(&key);
	commit_keygen(&key, rand);
	/* The permutation elements live in an LNP commitment, which needs its own
	 * key: it carries one message row per slot of the is_bin sub-proof. */
	lnp_keyinit(&lkey);
	lnp_keygen(&lkey, rand);

	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short_crt(r[i][j]);
		}
		commit_sample_short(m[i]);
		commit_init(&com[i]);
		commit_doit(&com[i], m[i], &key, r[i]);
	}

	/* Prover shuffles messages (only a circular shift for simplicity), and
	 * commits to the matching permutation elements sigma_i = g(pi(i)). */
	for (int i = 0; i < MSGS; i++) {
		pi[i] = (i + 1) % MSGS;
		nmod_poly_set(_m[i], m[pi[i]]);
		int_to_monomial(sigma[i], pi[i]);
	}

	TEST_ONCE("shuffle proof is consistent") {
		TEST_ASSERT(run(com, m, _m, sigma, r, &key, &lkey, rand) == 1, end);
	} TEST_END;

	/* Mount the attack of Section 4.1: the output list is obtained from the
	 * input list by swapping the first CRT component of two messages. It is
	 * therefore *not* a permutation of the input over R_p, yet the product
	 * identity that the original proof checked is still satisfied. */
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_set(am[i], m[i]);
	}
	crt_swap(am[0], am[1], m[0], m[1]);

	TEST_ONCE("CRT-mixed list is not a permutation but passes Neff's product") {
		TEST_ASSERT(nmod_poly_equal(am[0], m[0]) == 0, end);
		TEST_ASSERT(nmod_poly_equal(am[0], m[1]) == 0, end);
		commit_sample_rand(rho, rand, DEGCRT);
		neff_product(p0, m, rho);
		neff_product(p1, am, rho);
		TEST_ASSERT(nmod_poly_equal(p0, p1) == 1, end);
	} TEST_END;

	/* The cheating prover claims the identity permutation. Lemma 5 now ties
	 * every message to its index, so the two CRT components disagree and the
	 * product no longer matches. */
	for (int i = 0; i < MSGS; i++) {
		int_to_monomial(sigma[i], i);
	}

	TEST_ONCE("shuffle proof rejects the CRT-mixing attack") {
		TEST_ASSERT(run(com, m, am, sigma, r, &key, &lkey, rand) == 0, end);
	} TEST_END;

	/* Second-order attack, targeting the missing is_bin(sigma_i) sub-proof:
	 * the cheating prover applies to the index encodings the very same CRT
	 * swap it applied to the messages. Then sigma_0 = g(1) and sigma_1 = g(0)
	 * in the first CRT component, while sigma_i = g(i) in the second, so the
	 * product of Lemma 5 balances in both components again. The resulting
	 * sigma_0, sigma_1 are CRT-mixed, hence not short, which is exactly what
	 * the norm check on the committed sigma_i rules out. */
	int_to_monomial(g0, 0);
	int_to_monomial(g1, 1);
	crt_swap(sigma[0], sigma[1], g0, g1);

	TEST_ONCE("CRT-mixed sigma is outside D") {
		TEST_ASSERT(nmod_poly_equal(sigma[0], g0) == 0, end);
		TEST_ASSERT(nmod_poly_equal(sigma[0], g1) == 0, end);
		/* Honest permutation elements are short, the CRT-mixed one is not.
		 * The test uses norm2_leq rather than commit_norm2_sqr because the
		 * squared norm of a CRT-mixed element overflows a 64-bit integer. */
		TEST_ASSERT(commit_norm2_leq(g0, (uint64_t) DEGREE) == 1, end);
		TEST_ASSERT(commit_norm2_leq(g1, (uint64_t) DEGREE) == 1, end);
		TEST_ASSERT(commit_norm2_leq(sigma[0], (uint64_t) DEGREE) == 0, end);
		TEST_ASSERT(commit_norm2_leq(sigma[1], (uint64_t) DEGREE) == 0, end);
	} TEST_END;

	TEST_ONCE("shuffle proof rejects the CRT-mixing attack on sigma") {
		TEST_ASSERT(run(com, m, am, sigma, r, &key, &lkey, rand) == 0, end);
	} TEST_END;

  end:
	commit_finish();

	nmod_poly_clear(rho);
	nmod_poly_clear(p0);
	nmod_poly_clear(p1);
	nmod_poly_clear(g0);
	nmod_poly_clear(g1);
	for (int i = 0; i < MSGS; i++) {
		commit_free(&com[i]);
		nmod_poly_clear(m[i]);
		nmod_poly_clear(_m[i]);
		nmod_poly_clear(am[i]);
		nmod_poly_clear(sigma[i]);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(r[i][j][k]);
			}
		}
	}
	commit_keyfree(&key);
	lnp_keyfree(&lkey);
}

static void bench(flint_rand_t rand) {
	commitkey_t key;
	lnpkey_t lkey;
	commit_t com[MSGS];
	int pi[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS], sigma[MSGS];
	nmod_poly_t alpha, gamma, beta;
	nmod_poly_t r[MSGS][WIDTH][2];
	nmod_poly_t y[WIDTH][2], _y[WIDTH][2];
	nmod_poly_t t[2], tp[2], _t[2], u[2];

	nmod_poly_init(alpha, MODP);
	nmod_poly_init(gamma, MODP);
	nmod_poly_init(beta, MODP);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
		nmod_poly_init(sigma[i], MODP);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(r[i][j][k], MODP);
			}
		}
	}

	/* Generate commitment key-> */
	commit_setup();
	commit_keyinit(&key);
	commit_keygen(&key, rand);
	/* The permutation elements live in an LNP commitment, which needs its own
	 * key: it carries one message row per slot of the is_bin sub-proof. */
	lnp_keyinit(&lkey);
	lnp_keygen(&lkey, rand);

	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short_crt(r[i][j]);
		}
		commit_sample_short(m[i]);
		commit_init(&com[i]);
		commit_doit(&com[i], m[i], &key, r[i]);
	}

	/* Prover shuffles messages (only a circular shift for simplicity). */
	for (int i = 0; i < MSGS; i++) {
		pi[i] = (i + 1) % MSGS;
		nmod_poly_set(_m[i], m[pi[i]]);
		int_to_monomial(sigma[i], pi[i]);
	}

	BENCH_BEGIN("shuffle-proof (N messages)") {
		BENCH_ADD(run(com, m, _m, sigma, r, &key, &lkey, rand));
	} BENCH_END;

	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(t[i], MODP);
		nmod_poly_init(tp[i], MODP);
		nmod_poly_init(_t[i], MODP);
		nmod_poly_init(u[i], MODP);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_init(y[i][j], MODP);
			nmod_poly_init(_y[i][j], MODP);
		}
	}
	commit_finish();

	nmod_poly_clear(alpha);
	nmod_poly_clear(gamma);
	nmod_poly_clear(beta);
	for (int i = 0; i < MSGS; i++) {
		commit_free(&com[i]);
		nmod_poly_clear(m[i]);
		nmod_poly_clear(_m[i]);
		nmod_poly_clear(sigma[i]);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(r[i][j][k]);
			}
		}
	}
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(t[i]);
		nmod_poly_clear(tp[i]);
		nmod_poly_clear(_t[i]);
		nmod_poly_clear(u[i]);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_clear(y[i][j]);
			nmod_poly_clear(_y[i][j]);
		}
	}

	commit_keyfree(&key);
	lnp_keyfree(&lkey);
}

/* Select which phases to run: "test", "bench", or neither for both. Keeping
 * the benchmarks out of a test run matters in practice, since they dominate
 * the runtime by two orders of magnitude. */
static int phase_selected(int argc, char *argv[], const char *phase) {
	return argc < 2 || strcmp(argv[1], phase) == 0;
}

int main(int argc, char *argv[]) {
	flint_rand_t rand;

	flint_rand_init(rand);

	if (phase_selected(argc, argv, "test")) {
		printf("\n** Tests for lattice-based shuffle proof:\n\n");
		test(rand);
	}

	if (phase_selected(argc, argv, "bench")) {
		printf("\n** Benchmarks for lattice-based shuffle proof:\n\n");
		bench(rand);
	}

	flint_rand_clear(rand);
}
