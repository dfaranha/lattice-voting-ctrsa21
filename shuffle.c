#include <math.h>
#include <stdlib.h>

#include "param.h"
#include "commit.h"
#include "test.h"
#include "bench.h"
#include "serial.h"
#include "assert.h"
#include "fastrandombytes.h"
#include "randombytes.h"
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
 * Sample a uniformly random double in [0, 1) with 53 bits of precision.
 *
 * Rejection sampling must not reuse the Fiat-Shamir stream, so this draws from
 * the operating system rather than from fastrandombytes.
 *
 * @return the sampled value.
 */
static double uniform_double(void) {
	uint64_t bits;

	randombytes((unsigned char *)&bits, sizeof(bits));
	/* Keep the top 53 bits, the precision of the mantissa of a double. */
	return (double)(bits >> 11) * 0x1.0p-53;
}

/* z and v are in coefficient representation: the test needs their centred
 * coefficients, and the prover keeps them that way so as not to reconstruct
 * them from CRT form on every attempt. */
int rej_sampling(nmod_poly_t z[], nmod_poly_t v[], uint64_t s2, int width) {
	double r, u, M = 1.75;
	int64_t dot, norm;
	int64_t c0, c1;
	int result;

	u = uniform_double();

	norm = dot = 0;
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < DEGREE; j++) {
			c0 = nmod_poly_get_coeff_ui(z[i], j);
			c1 = nmod_poly_get_coeff_ui(v[i], j);
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
static void challenge_from_hash(nmod_poly_t d[2], nmod_poly_t dc,
		const uint8_t hash[SHA256HashSize]) {
	/* A local stream: keying the global one with public data would expose
	 * the masks that the Gaussian sampler draws from it. */
	prg_t prg;
	uint32_t buf;

	nmod_poly_t c;

	prg_setseed(&prg, hash);
	/* A challenge from the space C of the paper: NONZERO coefficients of +-1
	 * and the rest zero, so that ||d||_1 = NONZERO and |C| = 2^NONZERO *
	 * C(DEGREE, NONZERO). This used to be the difference of two binary
	 * polynomials, an element of C - C with ||d||_1 up to 2 * NONZERO, which
	 * put the extracted opening sqrt(2) above the bound binding is estimated
	 * at. DEGREE is a power of two, so the reduction below is unbiased, and
	 * the sign comes from the top bit, which the reduction does not use. */
	nmod_poly_init(c, MODP);
	nmod_poly_fit_length(c, DEGREE);
	for (int j = 0; j < NONZERO; j++) {
		uint32_t pos;

		do {
			prg_bytes(&prg, (unsigned char *)&buf, sizeof(buf));
			pos = buf % DEGREE;
		} while (nmod_poly_get_coeff_ui(c, pos) != 0);
		nmod_poly_set_coeff_ui(c, pos, (buf >> 31) ? MODP - 1 : 1);
	}
	pcrt_poly_reduce(d[0], c, 0);
	pcrt_poly_reduce(d[1], c, 1);
	if (dc != NULL) {
		nmod_poly_set(dc, c);
	}
	nmod_poly_clear(c);
}

void lin_hash(nmod_poly_t d[2], nmod_poly_t dc, commitkey_t *key,
		commit_t x, commit_t p,
		commit_t y, nmod_poly_t alpha, nmod_poly_t gamma, nmod_poly_t beta,
		nmod_poly_t u[2], nmod_poly_t t[2], nmod_poly_t tp[2],
		nmod_poly_t _t[2], nmod_poly_t vs[2], uint8_t *digest) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];

	SHA256Reset(&sha);

	/* Hash public key, through the digest commit_keygen computed once. */
	SHA256Input(&sha, key->digest, SHA256HashSize);

	/* Hash alpha, gamma, beta from linear relation. */
	commit_hash_poly(&sha, alpha);
	commit_hash_poly(&sha, gamma);
	commit_hash_poly(&sha, beta);

	/* Hash [x], [p], [x'], t, t_p, t' in CRT representation. */
	for (int i = 0; i < NCRT; i++) {
		commit_hash_crt(&sha, x.c1[i]);
		commit_hash_crt(&sha, x.c2[i]);
		commit_hash_crt(&sha, p.c1[i]);
		commit_hash_crt(&sha, p.c2[i]);
		commit_hash_crt(&sha, y.c1[i]);
		commit_hash_crt(&sha, y.c2[i]);
		commit_hash_crt(&sha, u[i]);
		commit_hash_crt(&sha, t[i]);
		commit_hash_crt(&sha, tp[i]);
		commit_hash_crt(&sha, _t[i]);
		commit_hash_crt(&sha, vs[i]);
	}

	SHA256Result(&sha, hash);
	if (digest != NULL) {
		memcpy(digest, hash, SHA256HashSize);
	}
	challenge_from_hash(d, dc, hash);
}

static void lin_prover(uint8_t digest[SHA256HashSize],
		nmod_poly_t y[WIDTH][2], nmod_poly_t w[WIDTH][2],
		nmod_poly_t _y[WIDTH][2], nmod_poly_t ys[1][2], nmod_poly_t t[2],
		nmod_poly_t tp[2], nmod_poly_t _t[2], nmod_poly_t vs[2],
		nmod_poly_t u[2], commit_t x, commit_t p, commit_t _x,
		commitkey_t *key, nmod_poly_t alpha, nmod_poly_t gamma,
		nmod_poly_t beta, nmod_poly_t r[WIDTH][2], nmod_poly_t s[WIDTH][2],
		nmod_poly_t _r[WIDTH][2], nmod_poly_t sig[2]) {
	/* The rejection tests work on coefficients, so the prover keeps the masks,
	 * the challenge and its products with the witness in coefficient form,
	 * and reduces into CRT form only what the commitments and the proof need.
	 * The witness is reconstructed once per proof rather than once per
	 * attempt. */
	nmod_poly_t tmp, d[2], dc, a[2], g[2];
	nmod_poly_t yc[WIDTH], wc[WIDTH], _yc[WIDTH], ysc[1];
	nmod_poly_t rc[WIDTH], sc[WIDTH], _rc[WIDTH], sigc;
	nmod_poly_t dr[WIDTH], ds[WIDTH], _dr[WIDTH], dsig[1];
	int rej0, rej1, rej2, rej3;
	/* Rejection sampling only makes z independent of the witness when it
	 * tests against the width y was actually sampled at, which is the SIGMA_C
	 * the Gaussian sampler is compiled with, not the paper's lower bound
	 * 11 * v * beta * sqrt(k * N) on that width. */
	const uint64_t sigma_sqr = (uint64_t) SIGMA_C * SIGMA_C;
	/* Masking the committed sigma is done with the much narrower SIGMA_S. */
	uint64_t sigma_s_sqr = (uint64_t) SIGMA_S * SIGMA_S;

	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(dr[i], MODP);
		nmod_poly_init(ds[i], MODP);
		nmod_poly_init(_dr[i], MODP);
		nmod_poly_init(yc[i], MODP);
		nmod_poly_init(wc[i], MODP);
		nmod_poly_init(_yc[i], MODP);
		nmod_poly_init(rc[i], MODP);
		nmod_poly_init(sc[i], MODP);
		nmod_poly_init(_rc[i], MODP);
		pcrt_poly_rec(rc[i], r[i]);
		pcrt_poly_rec(sc[i], s[i]);
		pcrt_poly_rec(_rc[i], _r[i]);
	}
	nmod_poly_init(dsig[0], MODP);
	nmod_poly_init(ysc[0], MODP);
	nmod_poly_init(sigc, MODP);
	pcrt_poly_rec(sigc, sig);
	nmod_poly_init(tmp, MODP);
	nmod_poly_init(dc, MODP);
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
			commit_sample_gauss(yc[i]);
			commit_sample_gauss(wc[i]);
			commit_sample_gauss(_yc[i]);
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_reduce(y[i][k], yc[i], k);
				pcrt_poly_reduce(w[i][k], wc[i], k);
				pcrt_poly_reduce(_y[i][k], _yc[i], k);
			}
		}
		/* Mask for the committed message sigma of p, narrow enough that the
		 * norm the verifier extracts from z_sigma stays below MODP^(1/2). */
		commit_sample_gauss_small(ysc[0]);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_reduce(ys[0][k], ysc[0], k);
		}
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
		lin_hash(d, dc, key, x, p, _x, alpha, gamma, beta, u, t, tp, _t, vs,
				digest);

		/* Prover */
		for (int i = 0; i < WIDTH; i++) {
			commit_poly_mulmod(dr[i], dc, rc[i]);
			nmod_poly_add(yc[i], yc[i], dr[i]);
			commit_poly_mulmod(ds[i], dc, sc[i]);
			nmod_poly_add(wc[i], wc[i], ds[i]);
			commit_poly_mulmod(_dr[i], dc, _rc[i]);
			nmod_poly_add(_yc[i], _yc[i], _dr[i]);
		}
		commit_poly_mulmod(dsig[0], dc, sigc);
		nmod_poly_add(ysc[0], ysc[0], dsig[0]);
		/* Each test is only worth running if all the previous ones accept. */
		rej0 = rej_sampling(yc, dr, sigma_sqr, WIDTH);
		rej1 = rej0 ? 1 : rej_sampling(wc, ds, sigma_sqr, WIDTH);
		rej2 = (rej0 || rej1) ? 1 : rej_sampling(_yc, _dr, sigma_sqr, WIDTH);
		if (rej0 || rej1 || rej2) {
			rej3 = 1;
		} else if (commit_norm2_leq(dsig[0],
						(uint64_t) DEGREE * SIGMA_S * SIGMA_S)) {
			rej3 = rej_sampling(ysc, dsig, sigma_s_sqr, 1);
		} else {
			/* If sigma lies outside D, d * sigma is far too large for the
			 * mask to hide and retrying cannot help. Stop rejecting and emit
			 * the transcript: the verifier rejects it on the norm check. */
			rej3 = 0;
		}
	} while (rej0 || rej1 || rej2 || rej3);

	/* The responses z = y + d * r, in the CRT form the proof carries. */
	for (int i = 0; i < WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_reduce(y[i][k], yc[i], k);
			pcrt_poly_reduce(w[i][k], wc[i], k);
			pcrt_poly_reduce(_y[i][k], _yc[i], k);
		}
		nmod_poly_clear(dr[i]);
		nmod_poly_clear(ds[i]);
		nmod_poly_clear(_dr[i]);
		nmod_poly_clear(yc[i]);
		nmod_poly_clear(wc[i]);
		nmod_poly_clear(_yc[i]);
		nmod_poly_clear(rc[i]);
		nmod_poly_clear(sc[i]);
		nmod_poly_clear(_rc[i]);
	}
	for (int k = 0; k < NCRT; k++) {
		pcrt_poly_reduce(ys[0][k], ysc[0], k);
	}
	nmod_poly_clear(dsig[0]);
	nmod_poly_clear(ysc[0]);
	nmod_poly_clear(sigc);
	nmod_poly_clear(dc);
	nmod_poly_clear(tmp);
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
	challenge_from_hash(_d, NULL, digest);

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
		lin_hash(ignored, NULL, key, x, p, _x, alpha, gamma, beta, u, t, tp, _t, vs,
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

/**
 * Derive the challenges X1 = tau and X2 = rho of Lemma 5 from the statement and
 * the commitments to the permutation elements.
 *
 * The product identity of Lemma 5 only shows that _m is a permutation of the
 * committed messages if tau and rho are fixed after the statement and after
 * the sigma_i, and cannot be chosen by the prover, which is what a verifier
 * drawing them at random gives an interactive proof. Here nobody is left to
 * draw them, so they are a hash of the input commitments, the output list and
 * the commitments P_i to the sigma_i -- everything the prover has fixed before
 * them, and all of which the verifier holds: the first two are the statement,
 * the P_i travel in the proof.
 *
 * tau must be invertible: where it vanishes in a CRT component, sigma_i
 * drops out of the product there, leaving the plain product of (m_i - rho)
 * that the CRT-mixing attack satisfies. That is a public condition, so a
 * counter absorbed into the hash re-derives both challenges in that case,
 * which keeps them deterministic. The other factor that must be invertible,
 * b_i = _m_i + sigma_i * tau - rho, involves the secret sigma_i, so only the
 * prover can test it, and it does so by committing to the sigma_i afresh,
 * which re-derives both challenges as well.
 *
 * @param[out] tau			- the challenge X1.
 * @param[out] rho			- the challenge X2.
 * @param[in] com			- the commitments to the input messages.
 * @param[in] p				- the commitments to the permutation elements.
 * @param[in] _m			- the output list.
 */
static void shuffle_challenges(nmod_poly_t tau, nmod_poly_t rho,
		commit_t com[MSGS], commit_t p[MSGS], nmod_poly_t _m[MSGS]) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	prg_t prg;

	for (uint8_t ctr = 0;; ctr++) {
		SHA256Reset(&sha);
		SHA256Input(&sha, (const uint8_t *)"tau,rho", 7);
		SHA256Input(&sha, &ctr, 1);
		for (int i = 0; i < MSGS; i++) {
			commit_hash_poly(&sha, _m[i]);
			for (int j = 0; j < NCRT; j++) {
				commit_hash_crt(&sha, com[i].c1[j]);
				commit_hash_crt(&sha, com[i].c2[j]);
				commit_hash_crt(&sha, p[i].c1[j]);
				commit_hash_crt(&sha, p[i].c2[j]);
			}
		}
		SHA256Result(&sha, hash);

		/* Of degree below N/2, so that differences of distinct challenges are
		 * invertible. */
		prg_setseed(&prg, hash);
		commit_sample_rand(tau, &prg, DEGCRT);
		commit_sample_rand(rho, &prg, DEGCRT);
		if (is_invertible(tau)) {
			return;
		}
	}
}

void shuffle_hash(nmod_poly_t beta, commit_t c[MSGS], commit_t p[MSGS],
		commit_t d[MSGS], nmod_poly_t _m[MSGS], nmod_poly_t tau,
		nmod_poly_t rho) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	prg_t prg;

	SHA256Reset(&sha);

	for (int i = 0; i < MSGS; i++) {
		commit_hash_poly(&sha, _m[i]);
		for (int j = 0; j < NCRT; j++) {
			commit_hash_crt(&sha, c[i].c1[j]);
			commit_hash_crt(&sha, c[i].c2[j]);
			commit_hash_crt(&sha, p[i].c1[j]);
			commit_hash_crt(&sha, p[i].c2[j]);
			commit_hash_crt(&sha, d[i].c1[j]);
			commit_hash_crt(&sha, d[i].c2[j]);
		}
	}
	commit_hash_poly(&sha, tau);
	commit_hash_poly(&sha, rho);
	SHA256Result(&sha, hash);

	/* beta is public and must be reproducible by the verifier, so it is drawn
	 * from a stream keyed by the hash rather than from the global one.
	 * Challenges are sampled of degree below DEGCRT = N/2, so that the
	 * difference of two distinct challenges is non-zero in both CRT components
	 * and hence invertible. */
	prg_setseed(&prg, hash);
	commit_sample_rand(beta, &prg, DEGCRT);
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
		nmod_poly_t r[MSGS][WIDTH][2], commitkey_t *key,
		uint8_t digest[MSGS][SHA256HashSize]) {
	int flag;
	nmod_poly_t beta, alpha, gamma, pub, tau, rho, t0, t1;
	nmod_poly_t a[MSGS], b[MSGS], theta[MSGS], inv[MSGS], _r[MSGS][WIDTH][2];
	nmod_poly_t pr[MSGS][WIDTH][2], sig[MSGS][2];

	nmod_poly_init(tau, MODP);
	nmod_poly_init(rho, MODP);
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
		nmod_poly_init(inv[i], MODP);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(sig[i][k], MODP);
			/* CRT representation of the committed permutation element. */
			pcrt_poly_reduce(sig[i][k], sigma[i], k);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_init(_r[i][j][k], MODP);
				nmod_poly_init(pr[i][j][k], MODP);
			}
		}
	}

	do {
		/* Prover commits to the permutation elements sigma_i of Lemma 5.
		 * This must happen *before* the challenges tau and rho: the soundness
		 * argument needs the sigma_i to be fixed independently of them, which
		 * deriving the challenges from the P_i enforces. */
		for (int i = 0; i < MSGS; i++) {
			for (int j = 0; j < WIDTH; j++) {
				commit_sample_short_crt(pr[i][j]);
			}
			commit_doit(&p[i], sigma[i], key, pr[i]);
		}
		shuffle_challenges(tau, rho, com, p, _m);

		/* Build the factors of the product from Lemma 5:
		 *     a_i = m_i + g(i) * tau - rho,
		 *     b_i = _m_i + sigma_i * tau - rho,
		 * with sigma_i = g(pi(i)) committed in P[i]. Tying the index to every
		 * message is what forces the permutations of the two CRT components
		 * to coincide, and is what the original product a_i - rho lacked. */
		flag = 1;
		for (int i = 0; i < MSGS; i++) {
			int_to_monomial(t1, i);
			commit_poly_mulmod(t0, t1, tau);
			nmod_poly_add(a[i], m[i], t0);
			nmod_poly_sub(a[i], a[i], rho);

			commit_poly_mulmod(t0, sigma[i], tau);
			nmod_poly_add(b[i], _m[i], t0);
			nmod_poly_sub(b[i], b[i], rho);
			/* The product argument inverts the b_i. This fails only with
			 * negligible probability, and then the prover commits to the
			 * sigma_i afresh, which re-derives the challenges. */
			flag &= is_invertible(b[i]);
		}
	} while (flag == 0);

	/* Prover samples theta_i and computes commitments D_i. */
	commit_sample_rand(theta[0], NULL, DEGREE);
	commit_poly_mulmod(t0, theta[0], b[0]);
	for (int j = 0; j < WIDTH; j++) {
		commit_sample_short_crt(_r[0][j]);
	}
	commit_doit(&d[0], t0, key, _r[0]);
	for (int i = 1; i < MSGS - 1; i++) {
		commit_sample_rand(theta[i], NULL, DEGREE);
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
	/* Every s_i divides by b_i. Inverting each is an extended gcd at full
	 * degree; inverting their product once and peeling the inverses off it
	 * (Montgomery's trick) costs one of those and three products per message.
	 * The loop above has already made sure every b_i is invertible, so their
	 * product is too. inv[i] holds the inverse of b[i]. */
	for (int i = 0; i < MSGS - 1; i++) {
		if (i == 0) {
			nmod_poly_set(inv[0], b[0]);
		} else {
			commit_poly_mulmod(inv[i], inv[i - 1], b[i]);
		}
	}
	nmod_poly_invmod(t0, inv[MSGS - 2], *commit_poly());
	for (int i = MSGS - 2; i > 0; i--) {
		/* t0 is the inverse of b[0] * ... * b[i]. */
		commit_poly_mulmod(inv[i], t0, inv[i - 1]);
		commit_poly_mulmod(t0, t0, b[i]);
	}
	nmod_poly_set(inv[0], t0);

	commit_poly_mulmod(s[0], theta[0], b[0]);
	commit_poly_mulmod(t0, beta, a[0]);
	nmod_poly_sub(s[0], s[0], t0);
	commit_poly_mulmod(s[0], s[0], inv[0]);
	for (int i = 1; i < MSGS - 1; i++) {
		commit_poly_mulmod(s[i], theta[i - 1], a[i]);
		commit_poly_mulmod(t0, theta[i], b[i]);
		nmod_poly_add(s[i], s[i], t0);
		commit_poly_mulmod(t0, s[i - 1], a[i]);
		nmod_poly_sub(s[i], s[i], t0);
		commit_poly_mulmod(s[i], s[i], inv[i]);
	}

	for (int l = 0; l < MSGS; l++) {
		shuffle_coeffs(alpha, gamma, pub, l, s, _m, beta, tau, rho);
		lin_prover(digest[l], y[l], w[l], _y[l], ys[l], t[l], tp[l], _t[l],
				vs[l], u[l], com[l], p[l], d[l], key, alpha, gamma, pub,
				r[l], pr[l], _r[l], sig[l]);
	}

	nmod_poly_clear(tau);
	nmod_poly_clear(rho);
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
		nmod_poly_clear(inv[i]);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(sig[i][k]);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_clear(_r[i][j][k]);
				nmod_poly_clear(pr[i][j][k]);
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
		nmod_poly_t _m[MSGS], commitkey_t *key,
		uint8_t digest[MSGS][SHA256HashSize]) {
	int result = 1;
	nmod_poly_t beta, alpha, gamma, pub, tau, rho;

	nmod_poly_init(tau, MODP);
	nmod_poly_init(rho, MODP);
	nmod_poly_init(beta, MODP);
	nmod_poly_init(alpha, MODP);
	nmod_poly_init(gamma, MODP);
	nmod_poly_init(pub, MODP);

	/* Verifier derives the challenges from the statement and the received P_i,
	 * as the prover did. */
	shuffle_challenges(tau, rho, com, p, _m);
	shuffle_hash(beta, com, p, d, _m, tau, rho);
	/* Now verify each \Prod_LIN instance, one for each commitment. */
	for (int l = 0; l < MSGS; l++) {
		shuffle_coeffs(alpha, gamma, pub, l, s, _m, beta, tau, rho);
		result &=
				lin_verifier(y[l], w[l], _y[l], ys[l], t[l], tp[l], _t[l],
				vs[l], u[l], com[l], p[l], d[l], key, alpha, gamma, pub,
				digest[l]);
	}

	nmod_poly_clear(tau);
	nmod_poly_clear(rho);
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
		 * than walked directly; the cost is the same DEGREE coefficients.
		 * There are MSGS - 1 of them: the last linear proof uses beta in
		 * place of s, so s[MSGS - 1] is neither computed nor checked, and
		 * sending it cost a ring element of zeros that anyone could alter. */
		if (l < MSGS - 1) {
			if (bw != NULL) {
				pcrt_poly_reduce(crt[0], pf->s[l], 0);
				pcrt_poly_reduce(crt[1], pf->s[l], 1);
				serial_put_uniform(bw, crt);
			} else {
				serial_get_uniform(r, crt);
				pcrt_poly_rec(pf->s[l], crt);
			}
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
		commitkey_t *key) {
	int result = 1;
	commit_t d[MSGS], p[MSGS];
	nmod_poly_t s[MSGS], u[MSGS][2];
	nmod_poly_t y[MSGS][WIDTH][2], w[MSGS][WIDTH][2], _y[MSGS][WIDTH][2];
	nmod_poly_t ys[MSGS][1][2], vs[MSGS][2];
	nmod_poly_t t[MSGS][2], tp[MSGS][2], _t[MSGS][2];

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
			}
		}
	}

	static uint8_t digest[MSGS][SHA256HashSize];

	/* The prover commits to the sigma_i and derives tau and rho from those
	 * commitments itself; the verifier re-derives them from the P_i it
	 * receives, so neither side is handed challenges from outside. */
	shuffle_prover(y, w, _y, ys, t, tp, _t, vs, u, d, s, com, p, m, _m, sigma,
			r, key, digest);

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

		/* Decode into objects of its own, so that a field the reader skipped
		 * cannot pass the prover's value on to the verifier, and reject what
		 * the reader flags. */
		commit_t d2[MSGS], p2[MSGS];
		nmod_poly_t s2[MSGS], ys2[MSGS][1][2];
		nmod_poly_t y2[MSGS][WIDTH][2], w2[MSGS][WIDTH][2], _y2[MSGS][WIDTH][2];
		static uint8_t digest2[MSGS][SHA256HashSize];

		for (int i = 0; i < MSGS; i++) {
			commit_init(&d2[i]);
			commit_init(&p2[i]);
			nmod_poly_init(s2[i], MODP);
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(ys2[i][0][k], MODP);
				for (int j = 0; j < WIDTH; j++) {
					nmod_poly_init(y2[i][j][k], MODP);
					nmod_poly_init(w2[i][j][k], MODP);
					nmod_poly_init(_y2[i][j][k], MODP);
				}
			}
		}
		memset(digest2, 0, sizeof(digest2));
		proof_t pf2 = { y2, w2, _y2, ys2, d2, p2, s2, digest2 };

		serial_reader_init(&br, buf, sizeof(buf));
		proof_walk(&pf2, NULL, &br);
		result = !br.overflow && shuffle_verifier(y2, w2, _y2, ys2, t, tp, _t,
				vs, u, d2, s2, com, p2, _m, key, digest2);

		for (int i = 0; i < MSGS; i++) {
			commit_free(&d2[i]);
			commit_free(&p2[i]);
			nmod_poly_clear(s2[i]);
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(ys2[i][0][k]);
				for (int j = 0; j < WIDTH; j++) {
					nmod_poly_clear(y2[i][j][k]);
					nmod_poly_clear(w2[i][j][k]);
					nmod_poly_clear(_y2[i][j][k]);
				}
			}
		}
	} else {
		result = shuffle_verifier(y, w, _y, ys, t, tp, _t, vs, u, d, s, com, p,
				_m, key, digest);
	}

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

static void test(void) {
	commitkey_t key;
	commit_t com[MSGS], com0[MSGS];
	int pi[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS], am[MSGS], sigma[MSGS], r[MSGS][WIDTH][2];
	nmod_poly_t m0[MSGS], _m0[MSGS], sigma0[MSGS];
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
		nmod_poly_init(m0[i], MODP);
		nmod_poly_init(_m0[i], MODP);
		nmod_poly_init(sigma0[i], MODP);
		commit_init(&com0[i]);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(r[i][j][k], MODP);
			}
		}
	}

	/* Generate commitment key-> */
	commit_setup();
	commit_keyinit(&key);
	commit_keygen(&key);

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

	TEST_ONCE("challenges have NONZERO coefficients of +-1") {
		nmod_poly_t d[2], c;
		uint8_t h[SHA256HashSize];
		int ones = 0, ternary = 1;

		nmod_poly_init(d[0], MODP);
		nmod_poly_init(d[1], MODP);
		nmod_poly_init(c, MODP);
		for (int i = 0; i < SHA256HashSize; i++) {
			h[i] = i;
		}
		challenge_from_hash(d, NULL, h);
		pcrt_poly_rec(c, d);
		for (int i = 0; i < DEGREE; i++) {
			ulong x = nmod_poly_get_coeff_ui(c, i);
			ternary &= (x <= 1 || x == MODP - 1);
			ones += (x != 0);
		}
		nmod_poly_clear(d[0]);
		nmod_poly_clear(d[1]);
		nmod_poly_clear(c);
		TEST_ASSERT(ternary && ones == NONZERO, end);
	} TEST_END;

	TEST_ONCE("decoder rejects encodings no writer produces") {
		static uint8_t raw[1 << 14], ok[1 << 14];
		bitreader_t br;
		bitwriter_t bw;
		nmod_poly_t a[2];
		int bad_u, bad_g, bad_e, good;

		nmod_poly_init(a[0], MODP);
		nmod_poly_init(a[1], MODP);
		/* All ones: a uniform coefficient of 2^32 - 1, at or above MODP, and a
		 * Gaussian escape carrying INT64_MIN, whose negation used to be
		 * undefined. */
		memset(raw, 0xFF, sizeof(raw));
		serial_reader_init(&br, raw, sizeof(raw));
		serial_get_uniform(&br, a);
		bad_u = br.overflow;
		serial_reader_init(&br, raw, sizeof(raw));
		serial_get_gauss(&br, a, SIGMA_C);
		bad_g = br.overflow;
		/* An escape carrying zero, which the unary code would have written. */
		memset(raw, 0, sizeof(raw));
		raw[0] = raw[1] = raw[2] = 0xFF;
		serial_reader_init(&br, raw, sizeof(raw));
		serial_get_gauss(&br, a, SIGMA_C);
		bad_e = br.overflow;
		/* What a writer produces reads back without the flag. */
		nmod_poly_zero(a[0]);
		nmod_poly_zero(a[1]);
		serial_writer_init(&bw, ok, sizeof(ok));
		serial_put_gauss(&bw, a, SIGMA_C);
		serial_put_uniform(&bw, a);
		serial_reader_init(&br, ok, sizeof(ok));
		serial_get_gauss(&br, a, SIGMA_C);
		serial_get_uniform(&br, a);
		good = !br.overflow;
		nmod_poly_clear(a[0]);
		nmod_poly_clear(a[1]);
		TEST_ASSERT(bad_u && bad_g && bad_e && good, end);
	} TEST_END;

	TEST_ONCE("proof survives serialisation, and measures what it should") {
		measure_proof = 1;
		TEST_ASSERT(run(com, m, _m, sigma, r, &key) == 1, end);
		measure_proof = 0;
		printf("\n    %zu bytes, %.1f KB per message\n", proof_bytes,
				proof_bytes / 1024.0 / MSGS);
		{
			/* The uniform part is a fixed width, but the Gaussian part is not
			 * one any more: it is coded against its distribution, so what it
			 * costs depends on the sample. The walk is therefore bracketed
			 * rather than predicted -- below what a flat encoding at the bound
			 * would have taken, above the entropy of the distribution it codes
			 * -- and a coding that stopped working would leave one side. Two
			 * widths here, since the opening of the permutation element is
			 * sampled at SIGMA_S, so the figure printed is their average. */
			double ncoef = 3.0 * WIDTH * DEGREE + DEGREE;
			/* Per message: D_i, P_i, and the digest; and s_i for all but the
			 * last. */
			double fixed = (4.0 * DEGREE * serial_bits_uniform()
					+ SHA256HashSize * 8)
					+ (MSGS - 1.0) / MSGS * DEGREE * serial_bits_uniform();
			double flat = fixed
					+ 3.0 * WIDTH * DEGREE * serial_bits_gauss(SIGMA_C)
					+ DEGREE * serial_bits_gauss(SIGMA_S);
			double ent = fixed
					+ 3.0 * WIDTH * DEGREE
					* log2(SIGMA_C * sqrt(2.0 * M_PI * M_E))
					+ DEGREE * log2(SIGMA_S * sqrt(2.0 * M_PI * M_E));

			TEST_ASSERT(proof_bytes * 8.0 < MSGS * flat, end);
			TEST_ASSERT(proof_bytes * 8.0 > MSGS * ent, end);
			printf("    %.2f bits a Gaussian coefficient, against %.2f flat "
					"and %.2f of entropy\n",
					(proof_bytes * 8.0 / MSGS - fixed) / ncoef,
					(flat - fixed) / ncoef, (ent - fixed) / ncoef);
		}
	} TEST_END;

	TEST_ONCE("shuffle proof is consistent") {
		TEST_ASSERT(run(com, m, _m, sigma, r, &key) == 1, end);
	} TEST_END;

	/* The proof must work on the caller's statement without altering it, so
	 * that a proof leaves nothing behind for the next one. */
	TEST_ONCE("shuffle proof leaves its statement untouched") {
		for (int i = 0; i < MSGS; i++) {
			nmod_poly_set(m0[i], m[i]);
			nmod_poly_set(_m0[i], _m[i]);
			nmod_poly_set(sigma0[i], sigma[i]);
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_set(com0[i].c1[k], com[i].c1[k]);
				nmod_poly_set(com0[i].c2[k], com[i].c2[k]);
			}
		}
		TEST_ASSERT(run(com, m, _m, sigma, r, &key) == 1, end);
		for (int i = 0; i < MSGS; i++) {
			TEST_ASSERT(nmod_poly_equal(m0[i], m[i]) == 1, end);
			TEST_ASSERT(nmod_poly_equal(_m0[i], _m[i]) == 1, end);
			TEST_ASSERT(nmod_poly_equal(sigma0[i], sigma[i]) == 1, end);
			for (int k = 0; k < NCRT; k++) {
				TEST_ASSERT(nmod_poly_equal(com0[i].c1[k], com[i].c1[k]), end);
				TEST_ASSERT(nmod_poly_equal(com0[i].c2[k], com[i].c2[k]), end);
			}
		}
	} TEST_END;

	TEST_ONCE("shuffle proof rejects a list that is not a permutation") {
		for (int i = 0; i < MSGS; i++) {
			nmod_poly_set(_m0[i], _m[i]);
		}
		/* Replace one output by a message that is not among the inputs, and
		 * keep the honest permutation elements. */
		nmod_poly_set(_m0[0], m[1]);
		nmod_poly_set_coeff_ui(_m0[0], 0,
				(nmod_poly_get_coeff_ui(_m0[0], 0) + 2) % MODP);
		TEST_ASSERT(run(com, m, _m0, sigma, r, &key) == 0, end);
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
		commit_sample_rand(rho, NULL, DEGCRT);
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
		TEST_ASSERT(run(com, m, am, sigma, r, &key) == 0, end);
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
		TEST_ASSERT(run(com, m, am, sigma, r, &key) == 0, end);
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
		commit_free(&com0[i]);
		nmod_poly_clear(m[i]);
		nmod_poly_clear(_m[i]);
		nmod_poly_clear(am[i]);
		nmod_poly_clear(sigma[i]);
		nmod_poly_clear(m0[i]);
		nmod_poly_clear(_m0[i]);
		nmod_poly_clear(sigma0[i]);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(r[i][j][k]);
			}
		}
	}
	commit_keyfree(&key);
}

static void bench(void) {
	commitkey_t key;
	commit_t com[MSGS];
	int pi[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS], sigma[MSGS];
	commit_t p, _com;
	nmod_poly_t alpha, gamma, beta, _mm, t0;
	nmod_poly_t r[MSGS][WIDTH][2], pr[WIDTH][2], _r[WIDTH][2];
	nmod_poly_t y[WIDTH][2], w[WIDTH][2], _y[WIDTH][2], ys[1][2];
	nmod_poly_t t[2], tp[2], _t[2], vs[2], u[2], sig[2];

	nmod_poly_init(alpha, MODP);
	nmod_poly_init(gamma, MODP);
	nmod_poly_init(beta, MODP);
	nmod_poly_init(_mm, MODP);
	nmod_poly_init(t0, MODP);
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
	commit_keygen(&key);

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
		BENCH_ADD(run(com, m, _m, sigma, r, &key));
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
			nmod_poly_init(pr[i][j], MODP);
			nmod_poly_init(_r[i][j], MODP);
		}
	}

	/* A statement the witness satisfies: com[0] opens to m[0] under r[0], p
	 * to the permutation element sigma[1] under randomness of its own, and
	 * _com to alpha * m[0] + gamma * sigma[1] + beta under randomness of its
	 * own. This used to pair com[0], com[1] and com[2] with freshly resampled
	 * randomness, and com[1] with sigma[1] it does not commit to: a false
	 * statement whose proof the verifier always rejected. */
	commit_sample_rand(alpha, NULL, DEGCRT);
	commit_sample_rand(gamma, NULL, DEGCRT);
	commit_sample_rand(beta, NULL, DEGREE);
	for (int i = 0; i < WIDTH; i++) {
		commit_sample_short_crt(pr[i]);
		commit_sample_short_crt(_r[i]);
	}
	commit_init(&p);
	commit_doit(&p, sigma[1], &key, pr);
	commit_poly_mulmod(_mm, alpha, m[0]);
	commit_poly_mulmod(t0, gamma, sigma[1]);
	nmod_poly_add(_mm, _mm, t0);
	nmod_poly_add(_mm, _mm, beta);
	commit_init(&_com);
	commit_doit(&_com, _mm, &key, _r);
	uint8_t dg[SHA256HashSize];

	BENCH_BEGIN("linear proof") {
		BENCH_ADD(lin_prover(dg, y, w, _y, ys, t, tp, _t, vs, u, com[0], p,
						_com, &key, alpha, gamma, beta, r[0], pr, _r, sig));
	} BENCH_END;

	BENCH_BEGIN("linear verifier") {
		BENCH_ADD(lin_verifier(y, w, _y, ys, t, tp, _t, vs, u, com[0], p,
						_com, &key, alpha, gamma, beta, dg));
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
			nmod_poly_clear(pr[i][j]);
			nmod_poly_clear(_r[i][j]);
		}
	}
	commit_free(&p);
	commit_free(&_com);
	nmod_poly_clear(_mm);
	nmod_poly_clear(t0);

	commit_keyfree(&key);
}

/* Select which phases to run: "test", "bench", or neither for both. Keeping
 * the benchmarks out of a test run matters in practice, since they dominate
 * the runtime by two orders of magnitude. */
static int phase_selected(int argc, char *argv[], const char *phase) {
	return argc < 2 || strcmp(argv[1], phase) == 0;
}

int main(int argc, char *argv[]) {
	if (phase_selected(argc, argv, "test")) {
		printf("\n** Tests for lattice-based shuffle proof:\n\n");
		test();
	}

	if (phase_selected(argc, argv, "bench")) {
		printf("\n** Benchmarks for lattice-based shuffle proof:\n\n");
		bench();
	}

	return test_failures() != 0;
}
