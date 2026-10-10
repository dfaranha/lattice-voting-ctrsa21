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
int rej_sampling(nmod_poly_t z[WIDTH], nmod_poly_t v[WIDTH], uint64_t s2) {
	double r, u, M = 1.75;
	int64_t dot, norm;
	int64_t c0, c1;
	int result;

	u = uniform_double();

	norm = dot = 0;
	for (int i = 0; i < WIDTH; i++) {
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

/* The challenge is a deterministic function of the digest, which is what lets
 * the proof carry the digest instead of the first messages: the verifier
 * derives the challenge from it, rebuilds each first message from the equation
 * that used to check it, and re-derives the digest to confirm it matches. */
static void challenge_from_hash(nmod_poly_t d[2], nmod_poly_t dc,
		const uint8_t hash[SHA256HashSize]) {
	uint32_t buf;
	/* A local stream, so that keying it with the public digest leaves the
	 * global one that produces the secret masks untouched. */
	prg_t prg;

	nmod_poly_t c;

	prg_setseed(&prg, (const uint8_t *) hash);
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

void lin_hash(nmod_poly_t d[2], nmod_poly_t dc, commitkey_t *key, commit_t x, commit_t y,
		nmod_poly_t alpha, nmod_poly_t beta, nmod_poly_t u[2],
		nmod_poly_t t[2], nmod_poly_t _t[2], uint8_t *digest) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];

	SHA256Reset(&sha);

	/* Hash public key, through the digest commit_keygen computed once. */
	SHA256Input(&sha, key->digest, SHA256HashSize);

	/* Hash alpha, beta from linear relation. */
	commit_hash_poly(&sha, alpha);
	commit_hash_poly(&sha, beta);

	/* Hash [x], [x'], t, t' in CRT representation. */
	for (int i = 0; i < NCRT; i++) {
		commit_hash_crt(&sha, x.c1[i]);
		commit_hash_crt(&sha, x.c2[i]);
		commit_hash_crt(&sha, y.c1[i]);
		commit_hash_crt(&sha, y.c2[i]);
		commit_hash_crt(&sha, u[i]);
		commit_hash_crt(&sha, t[i]);
		commit_hash_crt(&sha, _t[i]);
	}

	SHA256Result(&sha, hash);
	if (digest != NULL) {
		memcpy(digest, hash, SHA256HashSize);
	}
	challenge_from_hash(d, dc, hash);
}

static void lin_prover(uint8_t digest[SHA256HashSize],
		nmod_poly_t y[WIDTH][2], nmod_poly_t _y[WIDTH][2],
		nmod_poly_t t[2], nmod_poly_t _t[2], nmod_poly_t u[2],
		commit_t x, commit_t _x, commitkey_t *key, nmod_poly_t alpha,
		nmod_poly_t beta, nmod_poly_t r[WIDTH][2], nmod_poly_t _r[WIDTH][2],
		int l) {
	/* The rejection test works on coefficients, so the prover keeps the masks,
	 * the challenge and its product with the witness in coefficient form, and
	 * reduces into CRT form only what the commitments need. The witness is
	 * reconstructed once per proof rather than once per attempt. */
	nmod_poly_t tmp, d[2], dc, yc[WIDTH], _yc[WIDTH], rc[WIDTH], _rc[WIDTH];
	nmod_poly_t dr[WIDTH], _dr[WIDTH];
	int rej0, rej1;
	/* Rejection sampling only makes z independent of the witness when it
	 * tests against the width y was actually sampled at, which is the SIGMA_C
	 * the Gaussian sampler is compiled with, not the paper's lower bound
	 * 11 * v * beta * sqrt(k * N) on that width. */
	const uint64_t sigma_sqr = (uint64_t) SIGMA_C * SIGMA_C;

	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(dr[i], MODP);
		nmod_poly_init(_dr[i], MODP);
		nmod_poly_init(yc[i], MODP);
		nmod_poly_init(_yc[i], MODP);
		nmod_poly_init(rc[i], MODP);
		nmod_poly_init(_rc[i], MODP);
		pcrt_poly_rec(rc[i], r[i]);
		pcrt_poly_rec(_rc[i], _r[i]);
	}
	nmod_poly_init(tmp, MODP);
	nmod_poly_init(dc, MODP);
	nmod_poly_init(d[0], MODP);
	nmod_poly_init(d[1], MODP);

	do {
		for (int i = 0; i < NCRT; i++) {
			nmod_poly_zero(t[i]);
			nmod_poly_zero(_t[i]);
			nmod_poly_zero(u[i]);
			nmod_poly_zero(d[i]);
		}

		for (int i = 0; i < WIDTH; i++) {
			commit_sample_gauss(yc[i]);
			commit_sample_gauss(_yc[i]);
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_reduce(y[i][k], yc[i], k);
				pcrt_poly_reduce(_y[i][k], _yc[i], k);
			}
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
		}

		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < NCRT; j++) {
				pcrt_poly_mulmod(tmp, key->b2[i][j], y[i][j], j);
				pcrt_poly_mulmod(tmp, tmp, alpha, j);
				nmod_poly_add(u[j], u[j], tmp);
				pcrt_poly_mulmod(tmp, key->b2[i][j], _y[i][j], j);
				nmod_poly_sub(u[j], u[j], tmp);
			}
		}

		/* Sample challenge. */
		lin_hash(d, dc, key, x, _x, alpha, beta, u, t, _t, digest);

		/* Prover */
		for (int i = 0; i < WIDTH; i++) {
			commit_poly_mulmod(dr[i], dc, rc[i]);
			nmod_poly_add(yc[i], yc[i], dr[i]);
			commit_poly_mulmod(_dr[i], dc, _rc[i]);
			nmod_poly_add(_yc[i], _yc[i], _dr[i]);
		}
		/* The second test is only worth running if the first accepts. */
		rej0 = rej_sampling(yc, dr, sigma_sqr);
		rej1 = rej0 ? 1 : rej_sampling(_yc, _dr, sigma_sqr);
	} while (rej0 || rej1);

	/* The responses z = y + d * r, in the CRT form the proof carries. */
	for (int i = 0; i < WIDTH; i++) {
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_reduce(y[i][k], yc[i], k);
			pcrt_poly_reduce(_y[i][k], _yc[i], k);
		}
		nmod_poly_clear(dr[i]);
		nmod_poly_clear(_dr[i]);
		nmod_poly_clear(yc[i]);
		nmod_poly_clear(_yc[i]);
		nmod_poly_clear(rc[i]);
		nmod_poly_clear(_rc[i]);
	}
	nmod_poly_clear(dc);
	nmod_poly_clear(tmp);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(d[i]);
	}
}

static int lin_verifier(nmod_poly_t y[WIDTH][2], nmod_poly_t _y[WIDTH][2],
		nmod_poly_t t[2], nmod_poly_t _t[2], nmod_poly_t u[2],
		commit_t com, commit_t x, commitkey_t *key,
		nmod_poly_t alpha, nmod_poly_t beta, int l,
		const uint8_t digest[SHA256HashSize]) {
	nmod_poly_t tmp, _d[2], v[2], _v[2], z[WIDTH], _z[WIDTH], lin[2];
	int result = 1;

	nmod_poly_init(tmp, MODP);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(z[i], MODP);
		nmod_poly_init(_z[i], MODP);
	}
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(_d[i], MODP);
		nmod_poly_init(v[i], MODP);
		nmod_poly_init(_v[i], MODP);
		nmod_poly_init(lin[i], MODP);
		nmod_poly_zero(v[i]);
		nmod_poly_zero(_v[i]);
	}

	/* The proof carries the digest, not the first messages. */
	challenge_from_hash(_d, NULL, digest);

	/* Verifier checks norm, reconstruct from CRT representation. */
	for (int i = 0; i < WIDTH; i++) {
		pcrt_poly_rec(z[i], y[i]);
		pcrt_poly_rec(_z[i], _y[i]);
		/* Soundness checks: these must make verification fail rather than
		 * abort, and must not be compiled out by NDEBUG. */
		result &= commit_norm2_leq(z[i],
				(uint64_t) 4 * DEGREE * SIGMA_C * SIGMA_C);
		result &= commit_norm2_leq(_z[i],
				(uint64_t) 4 * DEGREE * SIGMA_C * SIGMA_C);
	}
	/* Verifier computes B1z and B1z'. */
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				pcrt_poly_mulmod(tmp, key->B1[i][j][k], y[j][k], k);
				nmod_poly_add(v[k], v[k], tmp);
				pcrt_poly_mulmod(tmp, key->B1[i][j][k], _y[j][k], k);
				nmod_poly_add(_v[k], _v[k], tmp);
			}
		}
	}
	/* The first messages are not transmitted. Each is recovered from the
	 * equation that used to check it, t = B_1z - d c1, and the check happens
	 * once at the end when the digest is rebuilt over them. */
	for (int j = 0; j < NCRT; j++) {
		pcrt_poly_mulmod(tmp, _d[j], com.c1[j], j);
		nmod_poly_sub(t[j], v[j], tmp);
		pcrt_poly_mulmod(tmp, _d[j], x.c1[j], j);
		nmod_poly_sub(_t[j], _v[j], tmp);
	}

	if (l == 0) {
		for (int j = 0; j < NCRT; j++) {
			pcrt_poly_mulmod(lin[j], alpha, com.c2[j], j);
			nmod_poly_add(lin[j], lin[j], beta);
			nmod_poly_sub(lin[j], lin[j], x.c2[j]);
			pcrt_poly_mulmod(lin[j], lin[j], _d[j], j);
			pcrt_poly_reduce(lin[j], lin[j], j);
		}
	}
	if (l > 0 && l < MSGS - 1) {
		for (int j = 0; j < NCRT; j++) {
			pcrt_poly_mulmod(lin[j], alpha, com.c2[j], j);
			nmod_poly_add(lin[j], lin[j], beta);
			nmod_poly_sub(lin[j], lin[j], x.c2[j]);
			pcrt_poly_mulmod(lin[j], lin[j], _d[j], j);
		}
	}
	if (l == MSGS - 1) {
		for (int j = 0; j < NCRT; j++) {
			pcrt_poly_mulmod(lin[j], alpha, com.c2[j], j);
			if (MSGS & 1) {
				nmod_poly_sub(lin[j], lin[j], beta);
			} else {
				nmod_poly_add(lin[j], lin[j], beta);
			}
			nmod_poly_sub(lin[j], lin[j], x.c2[j]);
			pcrt_poly_mulmod(lin[j], lin[j], _d[j], j);
		}
	}

	nmod_poly_zero(v[0]);
	nmod_poly_zero(v[1]);
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			pcrt_poly_mulmod(tmp, key->b2[i][j], y[i][j], j);
			pcrt_poly_mulmod(tmp, alpha, tmp, j);
			nmod_poly_add(v[j], v[j], tmp);
			pcrt_poly_mulmod(tmp, key->b2[i][j], _y[i][j], j);
			nmod_poly_sub(v[j], v[j], tmp);
		}
	}
	/* u = <b2 side> - d * <c2 side>, again a definition rather than a check. */
	for (int j = 0; j < NCRT; j++) {
		nmod_poly_sub(u[j], v[j], lin[j]);
	}

	/* One comparison replaces the three first-message checks: rebuilding the
	 * digest over the recovered first messages reproduces it exactly when each
	 * of those equations held, and not otherwise. */
	{
		uint8_t rebuilt[SHA256HashSize];
		nmod_poly_t ignored[2];

		for (int j = 0; j < NCRT; j++) {
			nmod_poly_init(ignored[j], MODP);
		}
		lin_hash(ignored, NULL, key, com, x, alpha, beta, u, t, _t, rebuilt);
		result &= (memcmp(rebuilt, digest, SHA256HashSize) == 0);
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_clear(ignored[j]);
		}
	}

	nmod_poly_clear(tmp);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_clear(z[i]);
		nmod_poly_clear(_z[i]);
	}
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(_d[i]);
		nmod_poly_clear(v[i]);
		nmod_poly_clear(_v[i]);
		nmod_poly_clear(lin[i]);
	}
	return result;
}

/**
 * Derive the evaluation point rho from the statement.
 *
 * The products of (m_i - rho) and (_m_i - rho) agreeing only shows that _m is
 * a permutation of m if rho is fixed after the statement and cannot be chosen
 * by the prover, which is what a verifier picking it at random gives an
 * interactive proof. Here nobody is left to pick it, so it is a hash of the
 * input commitments and the output list, which both sides compute alike. rho
 * must also differ from every _m_i, or the prover would have to invert zero;
 * a counter absorbed into the hash re-derives it in that case, which keeps it
 * deterministic.
 *
 * @param[out] rho			- the evaluation point.
 * @param[in] com			- the commitments to the input messages.
 * @param[in] _m			- the output list.
 */
static void shuffle_rho(nmod_poly_t rho, commit_t com[MSGS],
		nmod_poly_t _m[MSGS]) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	prg_t prg;
	int flag;

	for (uint8_t ctr = 0;; ctr++) {
		SHA256Reset(&sha);
		SHA256Input(&sha, (const uint8_t *)"rho", 3);
		SHA256Input(&sha, &ctr, 1);
		for (int i = 0; i < MSGS; i++) {
			commit_hash_poly(&sha, _m[i]);
			for (int j = 0; j < NCRT; j++) {
				commit_hash_crt(&sha, com[i].c1[j]);
				commit_hash_crt(&sha, com[i].c2[j]);
			}
		}
		SHA256Result(&sha, hash);

		prg_setseed(&prg, hash);
		commit_sample_rand(rho, &prg, DEGREE);
		flag = 1;
		for (int i = 0; i < MSGS; i++) {
			if (nmod_poly_equal(rho, _m[i]) == 1) {
				flag = 0;
			}
		}
		if (flag) {
			return;
		}
	}
}

/**
 * Shift the statement by rho into fresh copies, leaving the caller's intact.
 *
 * @param[out] sc			- the commitments with rho subtracted from c2.
 * @param[out] sm			- the list with rho subtracted from each entry.
 * @param[in] com			- the commitments.
 * @param[in] m				- the list.
 * @param[in] rho			- the evaluation point.
 */
static void shuffle_shift(commit_t sc[MSGS], nmod_poly_t sm[MSGS],
		commit_t com[MSGS], nmod_poly_t m[MSGS], nmod_poly_t rho) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_sub(sm[i], m[i], rho);
	}
	if (sc != NULL) {
		for (int j = 0; j < NCRT; j++) {
			pcrt_poly_reduce(t, rho, j);
			for (int i = 0; i < MSGS; i++) {
				nmod_poly_set(sc[i].c1[j], com[i].c1[j]);
				nmod_poly_sub(sc[i].c2[j], com[i].c2[j], t);
			}
		}
	}
	nmod_poly_clear(t);
}

void shuffle_hash(nmod_poly_t beta, commit_t c[MSGS], commit_t d[MSGS],
		nmod_poly_t _m[MSGS], nmod_poly_t rho) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	prg_t prg;

	SHA256Reset(&sha);

	for (int i = 0; i < MSGS; i++) {
		commit_hash_poly(&sha, _m[i]);
		for (int j = 0; j < NCRT; j++) {
			commit_hash_crt(&sha, c[i].c1[j]);
			commit_hash_crt(&sha, c[i].c2[j]);
			commit_hash_crt(&sha, d[i].c1[j]);
			commit_hash_crt(&sha, d[i].c2[j]);
		}
	}
	commit_hash_poly(&sha, rho);
	SHA256Result(&sha, hash);

	/* beta is public and must be reproducible by the verifier, so it is drawn
	 * from a stream keyed by the hash rather than from the global one. */
	prg_setseed(&prg, hash);
	commit_sample_rand(beta, &prg, DEGREE);
}

static void shuffle_prover(nmod_poly_t y[MSGS][WIDTH][2],
		nmod_poly_t _y[MSGS][WIDTH][2], nmod_poly_t t[MSGS][2],
		nmod_poly_t _t[MSGS][2], nmod_poly_t u[MSGS][2], commit_t d[MSGS],
		nmod_poly_t s[MSGS], commit_t com_in[MSGS], nmod_poly_t m_in[MSGS],
		nmod_poly_t _m_in[MSGS], nmod_poly_t r[MSGS][WIDTH][2],
		commitkey_t *key, uint8_t digest[MSGS][SHA256HashSize]) {
	nmod_poly_t beta, rho, t0, t1, theta[MSGS], inv[MSGS], _r[MSGS][WIDTH][2];
	/* The statement shifted by rho. */
	commit_t com[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS];

	nmod_poly_init(t0, MODP);
	nmod_poly_init(t1, MODP);
	nmod_poly_init(beta, MODP);
	nmod_poly_init(rho, MODP);
	for (int i = 0; i < MSGS; i++) {
		commit_init(&com[i]);
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
		nmod_poly_init(theta[i], MODP);
		nmod_poly_init(inv[i], MODP);
		for (int k = 0; k < NCRT; k++) {
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_init(_r[i][j][k], MODP);
			}
		}
	}

	/* Prover derives rho and shifts the commitments and both lists by it. */
	shuffle_rho(rho, com_in, _m_in);
	shuffle_shift(com, m, com_in, m_in, rho);
	shuffle_shift(NULL, _m, NULL, _m_in, rho);

	/* Prover samples theta_i and computes commitments D_i. */
	commit_sample_rand(theta[0], NULL, DEGREE);
	commit_poly_mulmod(t0, theta[0], _m[0]);
	for (int j = 0; j < WIDTH; j++) {
		commit_sample_short_crt(_r[0][j]);
	}
	commit_doit(&d[0], t0, key, _r[0]);
	for (int i = 1; i < MSGS - 1; i++) {
		commit_sample_rand(theta[i], NULL, DEGREE);
		commit_poly_mulmod(t0, theta[i - 1], m[i]);
		commit_poly_mulmod(t1, theta[i], _m[i]);
		nmod_poly_add(t0, t0, t1);
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short_crt(_r[i][j]);
		}
		commit_doit(&d[i], t0, key, _r[i]);
	}
	commit_poly_mulmod(t0, theta[MSGS - 2], m[MSGS - 1]);
	for (int j = 0; j < WIDTH; j++) {
		commit_sample_short_crt(_r[MSGS - 1][j]);
	}
	commit_doit(&d[MSGS - 1], t0, key, _r[MSGS - 1]);

	shuffle_hash(beta, com, d, _m, rho);
	/* Every s_i divides by _m_i. Inverting each is an extended gcd at full
	 * degree, which cost an eighth of the prover; inverting their product once
	 * and peeling the inverses off it (Montgomery's trick) costs one of those
	 * and three products per message. inv[i] holds the inverse of _m[i]. */
	for (int i = 0; i < MSGS - 1; i++) {
		if (i == 0) {
			nmod_poly_set(inv[0], _m[0]);
		} else {
			commit_poly_mulmod(inv[i], inv[i - 1], _m[i]);
		}
	}
	nmod_poly_invmod(t0, inv[MSGS - 2], *commit_poly());
	for (int i = MSGS - 2; i > 0; i--) {
		/* t0 is the inverse of _m[0] * ... * _m[i]. */
		commit_poly_mulmod(inv[i], t0, inv[i - 1]);
		commit_poly_mulmod(t0, t0, _m[i]);
	}
	nmod_poly_set(inv[0], t0);

	commit_poly_mulmod(s[0], theta[0], _m[0]);
	commit_poly_mulmod(t0, beta, m[0]);
	nmod_poly_sub(s[0], s[0], t0);
	commit_poly_mulmod(s[0], s[0], inv[0]);
	for (int i = 1; i < MSGS - 1; i++) {
		commit_poly_mulmod(s[i], theta[i - 1], m[i]);
		commit_poly_mulmod(t0, theta[i], _m[i]);
		nmod_poly_add(s[i], s[i], t0);
		commit_poly_mulmod(t0, s[i - 1], m[i]);
		nmod_poly_sub(s[i], s[i], t0);
		commit_poly_mulmod(s[i], s[i], inv[i]);
	}

	for (int l = 0; l < MSGS; l++) {
		if (l < MSGS - 1) {
			commit_poly_mulmod(t0, s[l], _m[l]);
		} else {
			commit_poly_mulmod(t0, beta, _m[l]);
		}

		if (l == 0) {
			lin_prover(digest[l], y[l], _y[l], t[l], _t[l], u[l], com[l], d[l],
					key, beta, t0, r[l], _r[l], l);
		} else {
			lin_prover(digest[l], y[l], _y[l], t[l], _t[l], u[l], com[l], d[l],
					key, s[l - 1], t0, r[l], _r[l], l);
		}
	}

	nmod_poly_clear(t0);
	nmod_poly_clear(t1);
	nmod_poly_clear(beta);
	nmod_poly_clear(rho);
	for (int i = 0; i < MSGS; i++) {
		commit_free(&com[i]);
		nmod_poly_clear(m[i]);
		nmod_poly_clear(_m[i]);
		nmod_poly_clear(theta[i]);
		nmod_poly_clear(inv[i]);
		for (int k = 0; k < NCRT; k++) {
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_clear(_r[i][j][k]);
			}
		}
	}
}

static int shuffle_verifier(nmod_poly_t y[MSGS][WIDTH][2],
		nmod_poly_t _y[MSGS][WIDTH][2], nmod_poly_t t[MSGS][2],
		nmod_poly_t _t[MSGS][2], nmod_poly_t u[MSGS][2], commit_t d[MSGS],
		nmod_poly_t s[MSGS], commit_t com_in[MSGS], nmod_poly_t _m_in[MSGS],
		commitkey_t *key, uint8_t digest[MSGS][SHA256HashSize]) {
	int result = 1;
	nmod_poly_t beta, rho, t0;
	/* The statement shifted by rho. */
	commit_t com[MSGS];
	nmod_poly_t _m[MSGS];

	nmod_poly_init(t0, MODP);
	nmod_poly_init(beta, MODP);
	nmod_poly_init(rho, MODP);
	for (int i = 0; i < MSGS; i++) {
		commit_init(&com[i]);
		nmod_poly_init(_m[i], MODP);
	}

	/* Verifier derives rho from the statement alone, as the prover did. */
	shuffle_rho(rho, com_in, _m_in);
	shuffle_shift(com, _m, com_in, _m_in, rho);

	shuffle_hash(beta, com, d, _m, rho);
	/* Now verify each \Prod_LIN instance, one for each commitment. */
	for (int l = 0; l < MSGS; l++) {
		if (l < MSGS - 1) {
			commit_poly_mulmod(t0, s[l], _m[l]);
		} else {
			commit_poly_mulmod(t0, beta, _m[l]);
		}

		if (l == 0) {
			result &=
					lin_verifier(y[l], _y[l], t[l], _t[l], u[l], com[l], d[l],
					key, beta, t0, l, digest[l]);
		} else {
			result &=
					lin_verifier(y[l], _y[l], t[l], _t[l], u[l], com[l], d[l],
					key, s[l - 1], t0, l, digest[l]);
		}
	}

	nmod_poly_clear(t0);
	nmod_poly_clear(beta);
	nmod_poly_clear(rho);
	for (int i = 0; i < MSGS; i++) {
		commit_free(&com[i]);
		nmod_poly_clear(_m[i]);
	}
	return result;
}

/*
 * Serialising the proof, so that its size is measured rather than modelled.
 *
 * Pointers to every part the prover sends, gathered so that one walk can both
 * pack and unpack them. Exactly one of the writer and the reader is non-NULL,
 * which is what keeps the two directions from drifting apart.
 *
 * What is absent matters as much as what is present. The commitments com[] and
 * the shuffled list _m[] are the statement, and the key is a public parameter,
 * so none of them is counted. Neither are the first messages t, _t and u: the
 * verifier recovers those from the equations, and the proof carries the
 * 32-byte digest in their place.
 */
typedef struct _proof_t {
	nmod_poly_t (*y)[WIDTH][2];
	nmod_poly_t (*_y)[WIDTH][2];
	commit_t *d;
	nmod_poly_t *s;
	uint8_t (*digest)[SHA256HashSize];
} proof_t;

static void walk_uniform(bitwriter_t *w, bitreader_t *r, nmod_poly_t a[2]) {
	if (w != NULL) {
		serial_put_uniform(w, a);
	} else {
		serial_get_uniform(r, a);
	}
}

static void walk_gauss(bitwriter_t *w, bitreader_t *r, nmod_poly_t a[2],
		ulong sigma) {
	if (w != NULL) {
		serial_put_gauss(w, a, sigma);
	} else {
		serial_get_gauss(r, a, sigma);
	}
}

static void proof_walk(proof_t *pf, bitwriter_t *w, bitreader_t *r) {
	nmod_poly_t crt[2];

	nmod_poly_init(crt[0], MODP);
	nmod_poly_init(crt[1], MODP);
	for (int l = 0; l < MSGS; l++) {
		for (int i = 0; i < WIDTH; i++) {
			walk_gauss(w, r, pf->y[l][i], SIGMA_C);
			walk_gauss(w, r, pf->_y[l][i], SIGMA_C);
		}
		walk_uniform(w, r, pf->d[l].c1);
		walk_uniform(w, r, pf->d[l].c2);
		/* s is held in coefficient representation, so it is converted rather
		 * than walked directly; the cost is the same DEGREE coefficients.
		 * There are MSGS - 1 of them: the last linear proof uses beta in
		 * place of s, so s[MSGS - 1] is neither computed nor checked, and
		 * sending it cost a ring element of zeros that anyone could alter. */
		if (l < MSGS - 1) {
			if (w != NULL) {
				pcrt_poly_reduce(crt[0], pf->s[l], 0);
				pcrt_poly_reduce(crt[1], pf->s[l], 1);
				serial_put_uniform(w, crt);
			} else {
				serial_get_uniform(r, crt);
				pcrt_poly_rec(pf->s[l], crt);
			}
		}
		for (int i = 0; i < SHA256HashSize; i++) {
			if (w != NULL) {
				serial_put_byte(w, pf->digest[l][i]);
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
		nmod_poly_t r[MSGS][WIDTH][2], commitkey_t *key) {
	int result = 1;
	commit_t d[MSGS];
	nmod_poly_t s[MSGS], u[MSGS][2];
	nmod_poly_t y[MSGS][WIDTH][2], _y[MSGS][WIDTH][2], t[MSGS][2], _t[MSGS][2];

	for (int i = 0; i < MSGS; i++) {
		commit_init(&d[i]);
		nmod_poly_init(s[i], MODP);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(t[i][k], MODP);
			nmod_poly_init(_t[i][k], MODP);
			nmod_poly_init(u[i][k], MODP);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_init(y[i][j][k], MODP);
				nmod_poly_init(_y[i][j][k], MODP);
			}
		}
	}

	static uint8_t digest[MSGS][SHA256HashSize];

	shuffle_prover(y, _y, t, _t, u, d, s, com, m, _m, r, key, digest);

	/* Round-trip the proof through its serialisation before verifying it, so
	 * that the measured size is the size of something that actually verifies
	 * and not just a count of struct fields. */
	if (measure_proof) {
		static uint8_t buf[4 << 20];
		bitwriter_t bw;
		bitreader_t br;
		proof_t pf = { y, _y, d, s, digest };

		serial_writer_init(&bw, buf, sizeof(buf));
		proof_walk(&pf, &bw, NULL);
		proof_bytes = bw.overflow ? 0 : serial_bytes(&bw);

		/* Decode into objects of its own, so that a field the reader skipped
		 * cannot pass the prover's value on to the verifier, and reject what
		 * the reader flags. */
		commit_t d2[MSGS];
		nmod_poly_t s2[MSGS], y2[MSGS][WIDTH][2], _y2[MSGS][WIDTH][2];
		static uint8_t digest2[MSGS][SHA256HashSize];

		for (int i = 0; i < MSGS; i++) {
			commit_init(&d2[i]);
			nmod_poly_init(s2[i], MODP);
			for (int j = 0; j < WIDTH; j++) {
				for (int k = 0; k < NCRT; k++) {
					nmod_poly_init(y2[i][j][k], MODP);
					nmod_poly_init(_y2[i][j][k], MODP);
				}
			}
		}
		memset(digest2, 0, sizeof(digest2));
		proof_t pf2 = { y2, _y2, d2, s2, digest2 };

		serial_reader_init(&br, buf, proof_bytes);
		proof_walk(&pf2, NULL, &br);
		serial_reader_finish(&br);
		result = !br.overflow && shuffle_verifier(y2, _y2, t, _t, u, d2, s2,
				com, _m, key, digest2);

		for (int i = 0; i < MSGS; i++) {
			commit_free(&d2[i]);
			nmod_poly_clear(s2[i]);
			for (int j = 0; j < WIDTH; j++) {
				for (int k = 0; k < NCRT; k++) {
					nmod_poly_clear(y2[i][j][k]);
					nmod_poly_clear(_y2[i][j][k]);
				}
			}
		}
	} else {
		result = shuffle_verifier(y, _y, t, _t, u, d, s, com, _m, key, digest);
	}

	for (int i = 0; i < MSGS; i++) {
		commit_free(&d[i]);
		nmod_poly_clear(s[i]);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(t[i][k]);
			nmod_poly_clear(_t[i][k]);
			nmod_poly_clear(u[i][k]);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_clear(y[i][j][k]);
				nmod_poly_clear(_y[i][j][k]);
			}
		}
	}

	return result;
}

static void test(void) {
	commitkey_t key;
	commit_t com[MSGS], com0[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS], m0[MSGS], _m0[MSGS], r[MSGS][WIDTH][2];

	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
		nmod_poly_init(m0[i], MODP);
		nmod_poly_init(_m0[i], MODP);
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

	/* Prover shuffles messages (only a circular shift for simplicity). */
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_set(_m[i], m[(i + 1) % MSGS]);
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
		nmod_poly_t a[2], t;
		int bad_u, bad_g, bad_e, good, exact, bad_t, bad_p;
		size_t len;

		nmod_poly_init(a[0], MODP);
		nmod_poly_init(a[1], MODP);
		nmod_poly_init(t, MODP);
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
		/* The reader must also end where the writer did: a byte past the
		 * end, or a padding bit set in the last byte, would give one proof
		 * a second byte string. A single coefficient of 2^j, for the first j
		 * whose code is not a whole number of bytes, leaves padding to set. */
		for (int j = 0; j < 30; j++) {
			nmod_poly_zero(t);
			nmod_poly_set_coeff_ui(t, 0, (ulong) 1 << j);
			pcrt_poly_reduce(a[0], t, 0);
			pcrt_poly_reduce(a[1], t, 1);
			serial_writer_init(&bw, ok, sizeof(ok));
			serial_put_gauss(&bw, a, SIGMA_C);
			if (bw.bits % 8 != 0) {
				break;
			}
		}
		len = serial_bytes(&bw);
		serial_reader_init(&br, ok, len);
		serial_get_gauss(&br, a, SIGMA_C);
		exact = (bw.bits % 8 != 0) && serial_reader_finish(&br);
		serial_reader_init(&br, ok, len + 1);
		serial_get_gauss(&br, a, SIGMA_C);
		bad_t = !serial_reader_finish(&br);
		ok[len - 1] |= 0x80;
		serial_reader_init(&br, ok, len);
		serial_get_gauss(&br, a, SIGMA_C);
		bad_p = !serial_reader_finish(&br);
		nmod_poly_clear(a[0]);
		nmod_poly_clear(a[1]);
		nmod_poly_clear(t);
		TEST_ASSERT(bad_u && bad_g && bad_e && good, end);
		TEST_ASSERT(exact && bad_t && bad_p, end);
	} TEST_END;

	TEST_ONCE("proof survives serialisation, and measures what it should") {
		measure_proof = 1;
		TEST_ASSERT(run(com, m, _m, r, &key) == 1, end);
		measure_proof = 0;
		printf("\n    %zu bytes, %.1f KB per message\n", proof_bytes,
				proof_bytes / 1024.0 / MSGS);
		{
			/* The uniform part is a fixed width, but the Gaussian part is not
			 * one any more: it is coded against its distribution, so what it
			 * costs depends on the sample. The walk is therefore bracketed
			 * rather than predicted -- below what a flat encoding at the bound
			 * would have taken, above the entropy of the distribution it codes
			 * -- and a coding that stopped working would leave one side. */
			/* Per message: D_i, two ring elements, and the digest; and s_i for
			 * all but the last. */
			double fixed = (2.0 * DEGREE * serial_bits_uniform()
					+ SHA256HashSize * 8)
					+ (MSGS - 1.0) / MSGS * DEGREE * serial_bits_uniform();
			double flat = fixed
					+ 2.0 * WIDTH * DEGREE * serial_bits_gauss(SIGMA_C);
			double ent = fixed + 2.0 * WIDTH * DEGREE
					* log2(SIGMA_C * sqrt(2.0 * M_PI * M_E));

			TEST_ASSERT(proof_bytes * 8.0 < MSGS * flat, end);
			TEST_ASSERT(proof_bytes * 8.0 > MSGS * ent, end);
			printf("    %.2f bits a Gaussian coefficient, against %d flat "
					"and %.2f of entropy\n",
					(proof_bytes * 8.0 / MSGS - fixed)
					/ (2.0 * WIDTH * DEGREE),
					serial_bits_gauss(SIGMA_C),
					log2(SIGMA_C * sqrt(2.0 * M_PI * M_E)));
		}
	} TEST_END;

	TEST_ONCE("shuffle proof is consistent") {
		TEST_ASSERT(run(com, m, _m, r, &key) == 1, end);
	} TEST_END;

	/* Shifting by rho used to happen in place, so every proof left the
	 * statement shifted for the next one. */
	TEST_ONCE("shuffle proof leaves its statement untouched") {
		for (int i = 0; i < MSGS; i++) {
			nmod_poly_set(m0[i], m[i]);
			nmod_poly_set(_m0[i], _m[i]);
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_set(com0[i].c1[k], com[i].c1[k]);
				nmod_poly_set(com0[i].c2[k], com[i].c2[k]);
			}
		}
		TEST_ASSERT(run(com, m, _m, r, &key) == 1, end);
		for (int i = 0; i < MSGS; i++) {
			TEST_ASSERT(nmod_poly_equal(m0[i], m[i]) == 1, end);
			TEST_ASSERT(nmod_poly_equal(_m0[i], _m[i]) == 1, end);
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
		/* Replace one output by a message that is not among the inputs. */
		nmod_poly_set(_m0[0], m[1]);
		nmod_poly_set_coeff_ui(_m0[0], 0,
				(nmod_poly_get_coeff_ui(_m0[0], 0) + 2) % MODP);
		TEST_ASSERT(run(com, m, _m0, r, &key) == 0, end);
	} TEST_END;

  end:
	commit_finish();

	for (int i = 0; i < MSGS; i++) {
		commit_free(&com[i]);
		commit_free(&com0[i]);
		nmod_poly_clear(m[i]);
		nmod_poly_clear(_m[i]);
		nmod_poly_clear(m0[i]);
		nmod_poly_clear(_m0[i]);
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
	nmod_poly_t m[MSGS], _m[MSGS];
	commit_t _com;
	nmod_poly_t alpha, beta, _mm;
	nmod_poly_t r[MSGS][WIDTH][2], _r[WIDTH][2], y[WIDTH][2], _y[WIDTH][2];
	nmod_poly_t t[2], _t[2], u[2];

	nmod_poly_init(alpha, MODP);
	nmod_poly_init(beta, MODP);
	nmod_poly_init(_mm, MODP);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
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
		nmod_poly_set(_m[i], m[(i + 1) % MSGS]);
	}

	BENCH_BEGIN("shuffle-proof (N messages)") {
		BENCH_ADD(run(com, m, _m, r, &key));
	} BENCH_END;

	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(t[i], MODP);
		nmod_poly_init(_t[i], MODP);
		nmod_poly_init(u[i], MODP);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_init(y[i][j], MODP);
			nmod_poly_init(_y[i][j], MODP);
			nmod_poly_init(_r[i][j], MODP);
		}
	}

	/* A statement the witness satisfies: com[0] opens to m[0] under r[0], and
	 * _com to alpha * m[0] + beta under randomness of its own. This used to
	 * pair com[0] and com[1] with r[0] for both and alpha, beta at random, a
	 * false statement whose proof the verifier always rejected. */
	commit_sample_rand(beta, NULL, DEGREE);
	commit_sample_rand(alpha, NULL, DEGREE);
	commit_poly_mulmod(_mm, alpha, m[0]);
	nmod_poly_add(_mm, _mm, beta);
	for (int i = 0; i < WIDTH; i++) {
		commit_sample_short_crt(_r[i]);
	}
	commit_init(&_com);
	commit_doit(&_com, _mm, &key, _r);
	uint8_t dg[SHA256HashSize];

	BENCH_BEGIN("linear proof") {
		BENCH_ADD(lin_prover(dg, y, _y, t, _t, u, com[0], _com, &key, alpha,
						beta, r[0], _r, 0));
	} BENCH_END;

	BENCH_BEGIN("linear verifier") {
		BENCH_ADD(lin_verifier(y, _y, t, _t, u, com[0], _com, &key, alpha,
						beta, 0, dg));
	} BENCH_END;

	commit_finish();

	nmod_poly_clear(alpha);
	nmod_poly_clear(beta);
	for (int i = 0; i < MSGS; i++) {
		commit_free(&com[i]);
		nmod_poly_clear(m[i]);
		nmod_poly_clear(_m[i]);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(r[i][j][k]);
			}
		}
	}
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(t[i]);
		nmod_poly_clear(_t[i]);
		nmod_poly_clear(u[i]);
	}
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < NCRT; j++) {
			nmod_poly_clear(y[i][j]);
			nmod_poly_clear(_y[i][j]);
			nmod_poly_clear(_r[i][j]);
		}
	}
	commit_free(&_com);
	nmod_poly_clear(_mm);

	commit_keyfree(&key);
}

/* Select which phases to run: "test", "bench", or neither for both. Keeping
 * the benchmarks out of a test run matters in practice, since they dominate
 * the runtime by two orders of magnitude. */
static int phase_selected(int argc, char *argv[], const char *phase) {
	return argc < 2 || strcmp(argv[1], phase) == 0;
}

int main(int argc, char *argv[]) {
	/* Unbuffered, so that a binary that dies still shows which test it
	 * was in: a crash used to take the buffered output with it. */
	setvbuf(stdout, NULL, _IONBF, 0);
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
