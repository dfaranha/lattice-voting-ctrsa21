#include <math.h>
#include <stdlib.h>

#include "param.h"
#include "commit.h"
#include "test.h"
#include "bench.h"
#include "serial.h"
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

/**
 * Sample a uniformly random double in [0, 1) with 53 bits of precision.
 *
 * Rejection sampling must not reuse the Fiat-Shamir stream, so this draws from
 * the operating system rather than from fastrandombytes.
 *
 * @return the sampled value.
 */
static double uniform_double(void) {
	uint64_t bits;

	getrandom(&bits, sizeof(bits), 0);
	/* Keep the top 53 bits, the precision of the mantissa of a double. */
	return (double)(bits >> 11) * 0x1.0p-53;
}

int rej_sampling(nmod_poly_t z[][2], nmod_poly_t v[][2], uint64_t s2,
		int width) {
	double r, u, M = 1.75;
	int64_t dot, norm;
	int64_t c0, c1;
	nmod_poly_t t0, t1;
	int result;

	nmod_poly_init(t0, MODP);
	nmod_poly_init(t1, MODP);

	u = uniform_double();

	norm = dot = 0;
	for (int i = 0; i < width; i++) {
		pcrt_poly_rec(t0, z[i]);
		pcrt_poly_rec(t1, v[i]);
		for (int j = 0; j < DEGREE; j++) {
			c0 = nmod_poly_get_coeff_ui(t0, j);
			c1 = nmod_poly_get_coeff_ui(t1, j);
			if (c0 > MODP / 2)
				c0 -= MODP;
			if (c1 > MODP / 2)
				c1 -= MODP;
			dot += c0 * c1;
			norm += c1 * c1;
		}
	}

	r = -2.0 * dot + norm;
	r = r / (2.0 * s2);
	r = exp(r) / M;

	result = u > r;

	nmod_poly_clear(t0);
	nmod_poly_clear(t1);
	return result;
}

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
/* The challenge is a deterministic function of the digest, which is what lets
 * the proof carry the digest instead of the first messages: the verifier
 * derives the challenge from it, rebuilds each first message from the equation
 * that used to check it, and re-derives the digest to confirm it matches. */
static void challenge_from_hash(nmod_poly_t d[2],
		const uint8_t hash[SHA256HashSize]) {
	uint32_t buf;

	fastrandombytes_setseed((uint8_t *) hash);
	/* The two slots of d are reused as the halves of the challenge d[0] - d[1]
	 * before being reduced into CRT components below, so this 2 belongs to the
	 * challenge construction rather than to the CRT decomposition. */
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

void lin_hash(nmod_poly_t d[2], commitkey_t *key, commit_t x, commit_t p,
		commit_t y, nmod_poly_t alpha, nmod_poly_t gamma, nmod_poly_t beta,
		nmod_poly_t u[2], nmod_poly_t t[2], nmod_poly_t tp[2],
		nmod_poly_t _t[2], nmod_poly_t vs[2], uint8_t *digest) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];

	SHA256Reset(&sha);

	/* Hash public key. */
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				hash_poly(&sha, key->B1[i][j][k]);
				if (i == 0) {
					hash_poly(&sha, key->b2[j][k]);
				}
			}
		}
	}

	/* Hash alpha, gamma, beta from linear relation. */
	hash_poly(&sha, alpha);
	hash_poly(&sha, gamma);
	hash_poly(&sha, beta);

	/* Hash [x], [p], [x'], t, t_p, t' in CRT representation. */
	for (int i = 0; i < NCRT; i++) {
		hash_poly(&sha, x.c1[i]);
		hash_poly(&sha, x.c2[i]);
		hash_poly(&sha, p.c1[i]);
		hash_poly(&sha, p.c2[i]);
		hash_poly(&sha, y.c1[i]);
		hash_poly(&sha, y.c2[i]);
		hash_poly(&sha, u[i]);
		hash_poly(&sha, t[i]);
		hash_poly(&sha, tp[i]);
		hash_poly(&sha, _t[i]);
		hash_poly(&sha, vs[i]);
	}

	SHA256Result(&sha, hash);
	if (digest != NULL) {
		memcpy(digest, hash, SHA256HashSize);
	}
	challenge_from_hash(d, hash);
}

static void lin_prover(uint8_t digest[SHA256HashSize],
		nmod_poly_t y[WIDTH][2], nmod_poly_t w[WIDTH][2],
		nmod_poly_t _y[WIDTH][2], nmod_poly_t ys[1][2], nmod_poly_t t[2],
		nmod_poly_t tp[2], nmod_poly_t _t[2], nmod_poly_t vs[2],
		nmod_poly_t u[2], commit_t x, commit_t p, commit_t _x,
		commitkey_t *key, nmod_poly_t alpha, nmod_poly_t gamma,
		nmod_poly_t beta, nmod_poly_t r[WIDTH][2], nmod_poly_t s[WIDTH][2],
		nmod_poly_t _r[WIDTH][2], nmod_poly_t sig[2]) {
	nmod_poly_t tmp, rec, d[2], a[2], g[2];
	nmod_poly_t dr[WIDTH][2], ds[WIDTH][2], _dr[WIDTH][2], dsig[1][2];
	int rej0, rej1, rej2, rej3;
	// Compute sigma^2 = (11 * v * beta * sqrt(k * N))^2.
	uint64_t sigma_sqr = 11 * NONZERO * BETA;
	sigma_sqr *= sigma_sqr * DEGREE * WIDTH;
	/* Masking the committed sigma is done with the much narrower SIGMA_S. */
	uint64_t sigma_s_sqr = (uint64_t) SIGMA_S * SIGMA_S;

	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_init(dr[i][j], MODP);
			nmod_poly_init(ds[i][j], MODP);
			nmod_poly_init(_dr[i][j], MODP);
		}
	}
	for (int j = 0; j < NCRT; j++) {
		nmod_poly_init(dsig[0][j], MODP);
	}
	nmod_poly_init(tmp, MODP);
	nmod_poly_init(rec, MODP);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(d[i], MODP);
		nmod_poly_init(a[i], MODP);
		nmod_poly_init(g[i], MODP);
		/* Reduce the public coefficients into CRT representation once. */
		pcrt_poly_reduce(a[i], alpha, i);
		pcrt_poly_reduce(g[i], gamma, i);
	}

	do {
		for (int i = 0; i < NCRT; i++) {
			nmod_poly_zero(t[i]);
			nmod_poly_zero(tp[i]);
			nmod_poly_zero(_t[i]);
			nmod_poly_zero(vs[i]);
			nmod_poly_zero(u[i]);
			nmod_poly_zero(d[i]);
		}

		for (int i = 0; i < WIDTH; i++) {
			commit_sample_gauss_crt(y[i]);
			commit_sample_gauss_crt(w[i]);
			commit_sample_gauss_crt(_y[i]);
		}
		/* Mask for the committed message sigma of p, narrow enough that the
		 * norm the verifier extracts from z_sigma stays below MODP^(1/2). */
		commit_sample_gauss_small_crt(ys[0]);
		for (int i = 0; i < HEIGHT; i++) {
			for (int j = 0; j < WIDTH; j++) {
				for (int k = 0; k < NCRT; k++) {
					pcrt_poly_mulmod(tmp, key->B1[i][j][k], y[j][k], k);
					nmod_poly_add(t[k], t[k], tmp);
					pcrt_poly_mulmod(tmp, key->B1[i][j][k], w[j][k], k);
					nmod_poly_add(tp[k], tp[k], tmp);
					pcrt_poly_mulmod(tmp, key->B1[i][j][k], _y[j][k], k);
					nmod_poly_add(_t[k], _t[k], tmp);
				}
			}
		}

		/* v_sigma = <b2, w> + y_sigma. */
		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < NCRT; j++) {
				pcrt_poly_mulmod(tmp, key->b2[i][j], w[i][j], j);
				nmod_poly_add(vs[j], vs[j], tmp);
			}
		}
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_add(vs[j], vs[j], ys[0][j]);
		}

		/* u = alpha * <b2, y> + gamma * <b2, w> - <b2, y'>. */
		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < NCRT; j++) {
				pcrt_poly_mulmod(tmp, key->b2[i][j], y[i][j], j);
				pcrt_poly_mulmod(tmp, tmp, a[j], j);
				nmod_poly_add(u[j], u[j], tmp);
				pcrt_poly_mulmod(tmp, key->b2[i][j], w[i][j], j);
				pcrt_poly_mulmod(tmp, tmp, g[j], j);
				nmod_poly_add(u[j], u[j], tmp);
				pcrt_poly_mulmod(tmp, key->b2[i][j], _y[i][j], j);
				nmod_poly_sub(u[j], u[j], tmp);
			}
		}

		/* Sample challenge. */
		lin_hash(d, key, x, p, _x, alpha, gamma, beta, u, t, tp, _t, vs,
				digest);

		/* Prover */
		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < NCRT; j++) {
				pcrt_poly_mulmod(dr[i][j], d[j], r[i][j], j);
				nmod_poly_add(y[i][j], y[i][j], dr[i][j]);
				pcrt_poly_mulmod(ds[i][j], d[j], s[i][j], j);
				nmod_poly_add(w[i][j], w[i][j], ds[i][j]);
				pcrt_poly_mulmod(_dr[i][j], d[j], _r[i][j], j);
				nmod_poly_add(_y[i][j], _y[i][j], _dr[i][j]);
			}
		}
		for (int j = 0; j < NCRT; j++) {
			pcrt_poly_mulmod(dsig[0][j], d[j], sig[j], j);
			nmod_poly_add(ys[0][j], ys[0][j], dsig[0][j]);
		}
		rej0 = rej_sampling(y, dr, sigma_sqr, WIDTH);
		rej1 = rej_sampling(w, ds, sigma_sqr, WIDTH);
		rej2 = rej_sampling(_y, _dr, sigma_sqr, WIDTH);
		/* If sigma lies outside D, d * sigma is far too large for the mask to
		 * hide and retrying cannot help. Stop rejecting and emit the
		 * transcript: the verifier rejects it on the norm check below. */
		pcrt_poly_rec(rec, dsig[0]);
		if (commit_norm2_leq(rec, (uint64_t) DEGREE * SIGMA_S * SIGMA_S)) {
			rej3 = rej_sampling(ys, dsig, sigma_s_sqr, 1);
		} else {
			rej3 = 0;
		}
	} while (rej0 || rej1 || rej2 || rej3);

	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_clear(dr[i][j]);
			nmod_poly_clear(ds[i][j]);
			nmod_poly_clear(_dr[i][j]);
		}
	}
	for (int j = 0; j < NCRT; j++) {
		nmod_poly_clear(dsig[0][j]);
	}
	nmod_poly_clear(tmp);
	nmod_poly_clear(rec);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(d[i]);
		nmod_poly_clear(a[i]);
		nmod_poly_clear(g[i]);
	}
}

static int lin_verifier(nmod_poly_t y[WIDTH][2], nmod_poly_t w[WIDTH][2],
		nmod_poly_t _y[WIDTH][2], nmod_poly_t ys[1][2], nmod_poly_t t[2],
		nmod_poly_t tp[2], nmod_poly_t _t[2], nmod_poly_t vs[2],
		nmod_poly_t u[2], commit_t x, commit_t p, commit_t _x,
		commitkey_t *key, nmod_poly_t alpha, nmod_poly_t gamma,
		nmod_poly_t beta, const uint8_t digest[SHA256HashSize]) {
	nmod_poly_t tmp, zs, _d[2], a[2], g[2], b[2];
	nmod_poly_t v[2], vp[2], _v[2], lhs[2], rhs[2];
	nmod_poly_t z[WIDTH], zp[WIDTH], _z[WIDTH];
	int result = 1;

	nmod_poly_init(tmp, MODP);
	nmod_poly_init(zs, MODP);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(z[i], MODP);
		nmod_poly_init(zp[i], MODP);
		nmod_poly_init(_z[i], MODP);
	}
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(_d[i], MODP);
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
	/* The proof carries the digest, not the first messages. Derive the
	 * challenge from it; the equations below rebuild each first message, and
	 * the digest is rebuilt over them at the end. */
	challenge_from_hash(_d, digest);

	/* Verifier checks norms, reconstructing from CRT representation. These are
	 * soundness checks and must make verification fail, so they are ordinary
	 * checks rather than assertions that vanish under NDEBUG. */
	for (int i = 0; i < WIDTH; i++) {
		pcrt_poly_rec(z[i], y[i]);
		pcrt_poly_rec(zp[i], w[i]);
		pcrt_poly_rec(_z[i], _y[i]);
		result &= commit_norm2_leq(z[i], (uint64_t) 4 * DEGREE * SIGMA_C * SIGMA_C);
		result &= commit_norm2_leq(zp[i], (uint64_t) 4 * DEGREE * SIGMA_C * SIGMA_C);
		result &= commit_norm2_leq(_z[i], (uint64_t) 4 * DEGREE * SIGMA_C * SIGMA_C);
	}

	/* The committed sigma must be short, which is what places it in the set D
	 * that Lemma 5 requires. This is the check that rules out a sigma mixing
	 * the two CRT components, which would otherwise let a cheating prover
	 * balance the product separately in each component. */
	pcrt_poly_rec(zs, ys[0]);
	result &= commit_norm2_leq(zs, (uint64_t) 4 * DEGREE * SIGMA_S * SIGMA_S);
	/* Verifier computes B1z, B1z_p and B1z'. */
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_mulmod(tmp, key->B1[i][j][k], y[j][k], k);
				nmod_poly_add(v[k], v[k], tmp);
				pcrt_poly_mulmod(tmp, key->B1[i][j][k], w[j][k], k);
				nmod_poly_add(vp[k], vp[k], tmp);
				pcrt_poly_mulmod(tmp, key->B1[i][j][k], _y[j][k], k);
				nmod_poly_add(_v[k], _v[k], tmp);
			}
		}
	}
	/* The first messages are not transmitted. Each is recovered from the
	 * equation that used to check it, t = B1z - d[x] and so on, and the check
	 * happens once at the end when the digest is rebuilt over them. */
	for (int j = 0; j < NCRT; j++) {
		pcrt_poly_mulmod(tmp, _d[j], x.c1[j], j);
		nmod_poly_sub(t[j], v[j], tmp);
		pcrt_poly_mulmod(tmp, _d[j], p.c1[j], j);
		nmod_poly_sub(tp[j], vp[j], tmp);
		pcrt_poly_mulmod(tmp, _d[j], _x.c1[j], j);
		nmod_poly_sub(_t[j], _v[j], tmp);
	}

	/* Verifier checks that <b2, z_p> + z_sigma = v_sigma + d * p.c2, which
	 * ties the short z_sigma to the message committed in p. */
	for (int j = 0; j < NCRT; j++) {
		pcrt_poly_mulmod(lhs[j], _d[j], p.c2[j], j);
		nmod_poly_zero(rhs[j]);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			pcrt_poly_mulmod(tmp, key->b2[i][j], w[i][j], j);
			nmod_poly_add(rhs[j], rhs[j], tmp);
		}
	}
	for (int j = 0; j < NCRT; j++) {
		nmod_poly_add(rhs[j], rhs[j], ys[0][j]);
		nmod_poly_sub(vs[j], rhs[j], lhs[j]);
	}

	/* Verifier checks the linear relation
	 *     alpha * <b2, z> + gamma * <b2, z_p> - <b2, z'>
	 *         = u + d * (alpha * x.c2 + gamma * p.c2 - x'.c2 + beta). */
	for (int j = 0; j < NCRT; j++) {
		pcrt_poly_mulmod(lhs[j], a[j], x.c2[j], j);
		pcrt_poly_mulmod(tmp, g[j], p.c2[j], j);
		nmod_poly_add(lhs[j], lhs[j], tmp);
		nmod_poly_sub(lhs[j], lhs[j], _x.c2[j]);
		nmod_poly_add(lhs[j], lhs[j], b[j]);
		pcrt_poly_mulmod(lhs[j], lhs[j], _d[j], j);
		nmod_poly_zero(rhs[j]);
	}

	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			pcrt_poly_mulmod(tmp, key->b2[i][j], y[i][j], j);
			pcrt_poly_mulmod(tmp, a[j], tmp, j);
			nmod_poly_add(rhs[j], rhs[j], tmp);
			pcrt_poly_mulmod(tmp, key->b2[i][j], w[i][j], j);
			pcrt_poly_mulmod(tmp, g[j], tmp, j);
			nmod_poly_add(rhs[j], rhs[j], tmp);
			pcrt_poly_mulmod(tmp, key->b2[i][j], _y[i][j], j);
			nmod_poly_sub(rhs[j], rhs[j], tmp);
		}
	}
	/* u = <b2 side> - d * <c2 side>, again a definition rather than a check. */
	for (int j = 0; j < NCRT; j++) {
		nmod_poly_sub(u[j], rhs[j], lhs[j]);
	}

	/* One comparison replaces the five first-message checks: rebuilding the
	 * digest over the recovered first messages reproduces it exactly when each
	 * of those equations held, and not otherwise. */
	{
		uint8_t rebuilt[SHA256HashSize];
		nmod_poly_t ignored[2];

		for (int j = 0; j < NCRT; j++) {
			nmod_poly_init(ignored[j], MODP);
		}
		lin_hash(ignored, key, x, p, _x, alpha, gamma, beta, u, t, tp, _t, vs,
				rebuilt);
		result &= (memcmp(rebuilt, digest, SHA256HashSize) == 0);
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_clear(ignored[j]);
		}
	}

	nmod_poly_clear(tmp);
	nmod_poly_clear(zs);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_clear(z[i]);
		nmod_poly_clear(zp[i]);
		nmod_poly_clear(_z[i]);
	}
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(_d[i]);
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

void shuffle_hash(nmod_poly_t beta, commit_t c[MSGS], commit_t p[MSGS],
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
			hash_poly(&sha, p[i].c1[j]);
			hash_poly(&sha, p[i].c2[j]);
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
		nmod_poly_t w[MSGS][WIDTH][2], nmod_poly_t _y[MSGS][WIDTH][2],
		nmod_poly_t ys[MSGS][1][2], nmod_poly_t t[MSGS][2],
		nmod_poly_t tp[MSGS][2], nmod_poly_t _t[MSGS][2],
		nmod_poly_t vs[MSGS][2], nmod_poly_t u[MSGS][2], commit_t d[MSGS],
		nmod_poly_t s[MSGS], commit_t com[MSGS], commit_t p[MSGS],
		nmod_poly_t m[MSGS], nmod_poly_t _m[MSGS], nmod_poly_t sigma[MSGS],
		nmod_poly_t r[MSGS][WIDTH][2], nmod_poly_t pr[MSGS][WIDTH][2],
		nmod_poly_t tau, nmod_poly_t rho, commitkey_t *key,
		uint8_t digest[MSGS][SHA256HashSize], flint_rand_t rng) {
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

	for (int l = 0; l < MSGS; l++) {
		shuffle_coeffs(alpha, gamma, pub, l, s, _m, beta, tau, rho);
		lin_prover(digest[l], y[l], w[l], _y[l], ys[l], t[l], tp[l], _t[l],
				vs[l], u[l], com[l], p[l], d[l], key, alpha, gamma, pub,
				r[l], pr[l], _r[l], sig[l]);
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
		nmod_poly_t w[MSGS][WIDTH][2], nmod_poly_t _y[MSGS][WIDTH][2],
		nmod_poly_t ys[MSGS][1][2], nmod_poly_t t[MSGS][2],
		nmod_poly_t tp[MSGS][2], nmod_poly_t _t[MSGS][2],
		nmod_poly_t vs[MSGS][2], nmod_poly_t u[MSGS][2], commit_t d[MSGS],
		nmod_poly_t s[MSGS], commit_t com[MSGS], commit_t p[MSGS],
		nmod_poly_t _m[MSGS], nmod_poly_t tau, nmod_poly_t rho,
		commitkey_t *key, uint8_t digest[MSGS][SHA256HashSize]) {
	int result = 1;
	nmod_poly_t beta, alpha, gamma, pub;

	nmod_poly_init(beta, MODP);
	nmod_poly_init(alpha, MODP);
	nmod_poly_init(gamma, MODP);
	nmod_poly_init(pub, MODP);

	shuffle_hash(beta, com, p, d, _m, tau, rho);
	/* Now verify each \Prod_LIN instance, one for each commitment. */
	for (int l = 0; l < MSGS; l++) {
		shuffle_coeffs(alpha, gamma, pub, l, s, _m, beta, tau, rho);
		result &=
				lin_verifier(y[l], w[l], _y[l], ys[l], t[l], tp[l], _t[l],
				vs[l], u[l], com[l], p[l], d[l], key, alpha, gamma, pub,
				digest[l]);
	}

	nmod_poly_clear(beta);
	nmod_poly_clear(alpha);
	nmod_poly_clear(gamma);
	nmod_poly_clear(pub);
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
/*
 * Serialising the proof, so that its size is measured rather than modelled.
 *
 * Pointers to every part the prover sends, gathered so that one walk can both
 * pack and unpack them. Exactly one of the writer and the reader is non-NULL,
 * which is what keeps the two directions from drifting apart.
 *
 * What is absent matters as much as what is present. The commitments com[] and
 * the shuffled list _m[] are the statement, and the key is a public parameter.
 * Neither are the first messages t, tp, _t, vs and u: the verifier recovers
 * those from the equations, and the proof carries the 32-byte digest instead.
 */
typedef struct _proof_t {
	nmod_poly_t (*y)[WIDTH][2];
	nmod_poly_t (*w)[WIDTH][2];
	nmod_poly_t (*_y)[WIDTH][2];
	nmod_poly_t (*ys)[1][2];
	commit_t *d;
	commit_t *p;
	nmod_poly_t *s;
	uint8_t (*digest)[SHA256HashSize];
} proof_t;

static void walk_uniform(bitwriter_t *bw, bitreader_t *r, nmod_poly_t a[2]) {
	if (bw != NULL) {
		serial_put_uniform(bw, a);
	} else {
		serial_get_uniform(r, a);
	}
}

static void walk_gauss(bitwriter_t *bw, bitreader_t *r, nmod_poly_t a[2],
		ulong sigma) {
	if (bw != NULL) {
		serial_put_gauss(bw, a, sigma);
	} else {
		serial_get_gauss(r, a, sigma);
	}
}

static void proof_walk(proof_t *pf, bitwriter_t *bw, bitreader_t *r) {
	nmod_poly_t crt[2];

	nmod_poly_init(crt[0], MODP);
	nmod_poly_init(crt[1], MODP);
	for (int l = 0; l < MSGS; l++) {
		for (int i = 0; i < WIDTH; i++) {
			walk_gauss(bw, r, pf->y[l][i], SIGMA_C);
			walk_gauss(bw, r, pf->w[l][i], SIGMA_C);
			walk_gauss(bw, r, pf->_y[l][i], SIGMA_C);
		}
		/* The short opening of the permutation element, at its own width. */
		walk_gauss(bw, r, pf->ys[l][0], SIGMA_S);
		walk_uniform(bw, r, pf->d[l].c1);
		walk_uniform(bw, r, pf->d[l].c2);
		walk_uniform(bw, r, pf->p[l].c1);
		walk_uniform(bw, r, pf->p[l].c2);
		/* s is held in coefficient representation, so it is converted rather
		 * than walked directly; the cost is the same DEGREE coefficients. */
		if (bw != NULL) {
			pcrt_poly_reduce(crt[0], pf->s[l], 0);
			pcrt_poly_reduce(crt[1], pf->s[l], 1);
			serial_put_uniform(bw, crt);
		} else {
			serial_get_uniform(r, crt);
			pcrt_poly_rec(pf->s[l], crt);
		}
		for (int i = 0; i < SHA256HashSize; i++) {
			if (bw != NULL) {
				serial_put_byte(bw, pf->digest[l][i]);
			} else {
				pf->digest[l][i] = serial_get_byte(r);
			}
		}
	}
	nmod_poly_clear(crt[0]);
	nmod_poly_clear(crt[1]);
}

/* Set by the size test: round-trip the proof through its serialisation before
 * verifying, and record how many bytes it took. */
static int measure_proof = 0;
static size_t proof_bytes = 0;

static int run(commit_t com[MSGS], nmod_poly_t m[MSGS], nmod_poly_t _m[MSGS],
		nmod_poly_t sigma[MSGS], nmod_poly_t r[MSGS][WIDTH][2],
		commitkey_t *key, flint_rand_t rng) {
	int result = 1;
	commit_t d[MSGS], p[MSGS];
	nmod_poly_t tau, rho, s[MSGS], u[MSGS][2];
	nmod_poly_t pr[MSGS][WIDTH][2];
	nmod_poly_t y[MSGS][WIDTH][2], w[MSGS][WIDTH][2], _y[MSGS][WIDTH][2];
	nmod_poly_t ys[MSGS][1][2], vs[MSGS][2];
	nmod_poly_t t[MSGS][2], tp[MSGS][2], _t[MSGS][2];

	nmod_poly_init(tau, MODP);
	nmod_poly_init(rho, MODP);
	for (int i = 0; i < MSGS; i++) {
		commit_init(&d[i]);
		commit_init(&p[i]);
		nmod_poly_init(s[i], MODP);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(t[i][k], MODP);
			nmod_poly_init(tp[i][k], MODP);
			nmod_poly_init(_t[i][k], MODP);
			nmod_poly_init(ys[i][0][k], MODP);
			nmod_poly_init(vs[i][k], MODP);
			nmod_poly_init(u[i][k], MODP);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_init(y[i][j][k], MODP);
				nmod_poly_init(w[i][j][k], MODP);
				nmod_poly_init(_y[i][j][k], MODP);
				nmod_poly_init(pr[i][j][k], MODP);
			}
		}
	}

	/* Prover commits to the permutation elements sigma_i of Lemma 5. This must
	 * happen *before* the challenges tau and rho are drawn: the soundness
	 * argument needs the sigma_i to be fixed independently of them. */
	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short_crt(pr[i][j]);
		}
		commit_doit(&p[i], sigma[i], key, pr[i]);
	}

	/* Verifier samples the two challenges X1 = tau and X2 = rho of Lemma 5,
	 * of degree below N/2 so that differences of distinct challenges are
	 * invertible. */
	commit_sample_rand(tau, rng, DEGCRT);
	commit_sample_rand(rho, rng, DEGCRT);

	static uint8_t digest[MSGS][SHA256HashSize];

	shuffle_prover(y, w, _y, ys, t, tp, _t, vs, u, d, s, com, p, m, _m, sigma,
			r, pr, tau, rho, key, digest, rng);

	/* Round-trip the proof through its serialisation before verifying it, so
	 * that the measured size is the size of something that actually verifies
	 * and not just a count of struct fields. */
	if (measure_proof) {
		static uint8_t buf[8 << 20];
		bitwriter_t bw;
		bitreader_t br;
		proof_t pf = { y, w, _y, ys, d, p, s, digest };

		serial_writer_init(&bw, buf, sizeof(buf));
		proof_walk(&pf, &bw, NULL);
		proof_bytes = bw.overflow ? 0 : serial_bytes(&bw);
		serial_reader_init(&br, buf, sizeof(buf));
		proof_walk(&pf, NULL, &br);
	}

	result = shuffle_verifier(y, w, _y, ys, t, tp, _t, vs, u, d, s, com, p, _m,
			tau, rho, key, digest);

	nmod_poly_clear(tau);
	nmod_poly_clear(rho);
	for (int i = 0; i < MSGS; i++) {
		commit_free(&d[i]);
		commit_free(&p[i]);
		nmod_poly_clear(s[i]);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(t[i][k]);
			nmod_poly_clear(tp[i][k]);
			nmod_poly_clear(_t[i][k]);
			nmod_poly_clear(ys[i][0][k]);
			nmod_poly_clear(vs[i][k]);
			nmod_poly_clear(u[i][k]);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_clear(y[i][j][k]);
				nmod_poly_clear(w[i][j][k]);
				nmod_poly_clear(_y[i][j][k]);
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

	TEST_ONCE("proof survives serialisation, and measures what it should") {
		measure_proof = 1;
		TEST_ASSERT(run(com, m, _m, sigma, r, &key, rand) == 1, end);
		measure_proof = 0;
		printf("\n    %zu bytes, %.1f KB per message\n", proof_bytes,
				proof_bytes / 1024.0 / MSGS);
		{
			size_t bu = serial_bits_uniform();
			size_t per = 5 * DEGREE * bu
					+ 3 * WIDTH * DEGREE * serial_bits_gauss(SIGMA_C)
					+ DEGREE * serial_bits_gauss(SIGMA_S)
					+ SHA256HashSize * 8;

			TEST_ASSERT(proof_bytes == (MSGS * per + 7) / 8, end);
		}
	} TEST_END;

	TEST_ONCE("shuffle proof is consistent") {
		TEST_ASSERT(run(com, m, _m, sigma, r, &key, rand) == 1, end);
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
		TEST_ASSERT(run(com, m, am, sigma, r, &key, rand) == 0, end);
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
		TEST_ASSERT(run(com, m, am, sigma, r, &key, rand) == 0, end);
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
}

static void bench(flint_rand_t rand) {
	commitkey_t key;
	commit_t com[MSGS];
	int pi[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS], sigma[MSGS];
	nmod_poly_t alpha, gamma, beta;
	nmod_poly_t r[MSGS][WIDTH][2];
	nmod_poly_t y[WIDTH][2], w[WIDTH][2], _y[WIDTH][2], ys[1][2];
	nmod_poly_t t[2], tp[2], _t[2], vs[2], u[2], sig[2];

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
		BENCH_ADD(run(com, m, _m, sigma, r, &key, rand));
	} BENCH_END;

	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(t[i], MODP);
		nmod_poly_init(tp[i], MODP);
		nmod_poly_init(_t[i], MODP);
		nmod_poly_init(vs[i], MODP);
		nmod_poly_init(u[i], MODP);
		nmod_poly_init(ys[0][i], MODP);
		nmod_poly_init(sig[i], MODP);
		pcrt_poly_reduce(sig[i], sigma[1], i);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_init(y[i][j], MODP);
			nmod_poly_init(w[i][j], MODP);
			nmod_poly_init(_y[i][j], MODP);
		}
	}

	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short_crt(r[i][j]);
		}
	}
	commit_sample_rand(alpha, rand, DEGCRT);
	commit_sample_rand(gamma, rand, DEGCRT);
	commit_sample_rand(beta, rand, DEGREE);
	uint8_t dg[SHA256HashSize];

	BENCH_BEGIN("linear proof") {
		BENCH_ADD(lin_prover(dg, y, w, _y, ys, t, tp, _t, vs, u, com[0],
						com[1], com[2], &key, alpha, gamma, beta, r[0], r[1],
						r[2], sig));
	} BENCH_END;

	BENCH_BEGIN("linear verifier") {
		BENCH_ADD(lin_verifier(y, w, _y, ys, t, tp, _t, vs, u, com[0], com[1],
						com[2], &key, alpha, gamma, beta, dg));
	} BENCH_END;

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
		nmod_poly_clear(vs[i]);
		nmod_poly_clear(u[i]);
		nmod_poly_clear(ys[0][i]);
		nmod_poly_clear(sig[i]);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_clear(y[i][j]);
			nmod_poly_clear(w[i][j]);
			nmod_poly_clear(_y[i][j]);
		}
	}

	commit_keyfree(&key);
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
