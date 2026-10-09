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

int rej_sampling(nmod_poly_t z[WIDTH], nmod_poly_t v[WIDTH], uint64_t s2) {
	double r, u, M = 1.75;
	int64_t dot, norm;
	int64_t c0, c1;
	int result;

	u = uniform_double();

	norm = dot = 0;
	/* No reconstruction: z and v are already in coefficient representation. */
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
static void challenge_from_hash(nmod_poly_t d,
		const uint8_t hash[SHA256HashSize]) {
	/* A local stream: keying the global one with public data would expose
	 * the masks that the Gaussian sampler draws from it. */
	prg_t prg;
	uint32_t buf;

	prg_setseed(&prg, hash);
	/* The challenge is the difference of two polynomials with NONZERO ones.
	 * The 2 here belongs to that construction, not to any CRT splitting. */
	nmod_poly_t c[2];
	for (int i = 0; i < 2; i++) {
		nmod_poly_init(c[i], MODP);
		nmod_poly_fit_length(c[i], DEGREE);
		for (int j = 0; j < NONZERO; j++) {
			prg_bytes(&prg, (unsigned char *)&buf, sizeof(buf));
			buf = buf % DEGREE;
			while (nmod_poly_get_coeff_ui(c[i], buf) != 0) {
				prg_bytes(&prg, (unsigned char *)&buf, sizeof(buf));
				buf = buf % DEGREE;
			}
			nmod_poly_set_coeff_ui(c[i], buf, 1);
		}
	}
	nmod_poly_sub(d, c[0], c[1]);
	nmod_poly_clear(c[0]);
	nmod_poly_clear(c[1]);
}

void lin_hash(nmod_poly_t d, commitkey_t *key, commit_t x, commit_t y,
		nmod_poly_t alpha, nmod_poly_t beta, nmod_poly_t u,
		nmod_poly_t t, nmod_poly_t _t, uint8_t *digest) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];

	SHA256Reset(&sha);

	/* Hash public key. */
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			SHA256Input(&sha, (const uint8_t *)key->B1[i][j]->coeffs,
					key->B1[i][j]->alloc * sizeof(uint64_t));
			if (i == 0) {
				SHA256Input(&sha, (const uint8_t *)key->b2[j]->coeffs,
						key->b2[j]->alloc * sizeof(uint64_t));
			}
		}
	}

	/* Hash alpha, beta from linear relation. */
	SHA256Input(&sha, (const uint8_t *)alpha->coeffs,
			alpha->alloc * sizeof(uint64_t));
	SHA256Input(&sha, (const uint8_t *)beta->coeffs,
			beta->alloc * sizeof(uint64_t));

	/* Hash [x], [x'], t, t'. */
	SHA256Input(&sha, (const uint8_t *)x.c1->coeffs,
			x.c1->alloc * sizeof(uint64_t));
	SHA256Input(&sha, (const uint8_t *)x.c2->coeffs,
			x.c2->alloc * sizeof(uint64_t));
	SHA256Input(&sha, (const uint8_t *)y.c1->coeffs,
			y.c1->alloc * sizeof(uint64_t));
	SHA256Input(&sha, (const uint8_t *)y.c2->coeffs,
			y.c2->alloc * sizeof(uint64_t));
	SHA256Input(&sha, (const uint8_t *)u->coeffs, u->alloc * sizeof(uint64_t));
	SHA256Input(&sha, (const uint8_t *)t->coeffs, t->alloc * sizeof(uint64_t));
	SHA256Input(&sha, (const uint8_t *)_t->coeffs,
			_t->alloc * sizeof(uint64_t));

	SHA256Result(&sha, hash);
	if (digest != NULL) {
		memcpy(digest, hash, SHA256HashSize);
	}
	challenge_from_hash(d, hash);
}

static void lin_prover(uint8_t digest[SHA256HashSize],
		nmod_poly_t y[WIDTH], nmod_poly_t _y[WIDTH],
		nmod_poly_t t, nmod_poly_t _t, nmod_poly_t u,
		commit_t x, commit_t _x, commitkey_t *key, nmod_poly_t alpha,
		nmod_poly_t beta, nmod_poly_t r[WIDTH], nmod_poly_t _r[WIDTH],
		int l) {
	nmod_poly_t tmp, d, dr[WIDTH], _dr[WIDTH];
	int rej0, rej1;
	/* Rejection sampling only makes z independent of the witness when it
	 * tests against the width y was actually sampled at, which is the SIGMA_C
	 * the Gaussian sampler is compiled with, not the paper's lower bound
	 * 11 * v * beta * sqrt(k * N) on that width. */
	const uint64_t sigma_sqr = (uint64_t) SIGMA_C * SIGMA_C;

	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(dr[i], MODP);
		nmod_poly_init(_dr[i], MODP);
	}
	nmod_poly_init(tmp, MODP);
	nmod_poly_init(d, MODP);

	do {
		nmod_poly_zero(t);
		nmod_poly_zero(_t);
		nmod_poly_zero(u);
		nmod_poly_zero(d);

		for (int i = 0; i < WIDTH; i++) {
			commit_sample_gauss(y[i]);
			commit_sample_gauss(_y[i]);
		}
		for (int i = 0; i < HEIGHT; i++) {
			for (int j = 0; j < WIDTH; j++) {
				commit_poly_mulmod(tmp, key->B1[i][j], y[j]);
				nmod_poly_add(t, t, tmp);
				commit_poly_mulmod(tmp, key->B1[i][j], _y[j]);
				nmod_poly_add(_t, _t, tmp);
			}
		}

		for (int i = 0; i < WIDTH; i++) {
			commit_poly_mulmod(tmp, key->b2[i], y[i]);
			commit_poly_mulmod(tmp, tmp, alpha);
			nmod_poly_add(u, u, tmp);
			commit_poly_mulmod(tmp, key->b2[i], _y[i]);
			nmod_poly_sub(u, u, tmp);
		}

		/* Sample challenge. */
		lin_hash(d, key, x, _x, alpha, beta, u, t, _t, digest);

		/* Prover */
		for (int i = 0; i < WIDTH; i++) {
			commit_poly_mulmod(dr[i], d, r[i]);
			nmod_poly_add(y[i], y[i], dr[i]);
			commit_poly_mulmod(_dr[i], d, _r[i]);
			nmod_poly_add(_y[i], _y[i], _dr[i]);
		}
		rej0 = rej_sampling(y, dr, sigma_sqr);
		rej1 = rej_sampling(_y, _dr, sigma_sqr);
	} while (rej0 || rej1);

	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_clear(dr[i]);
		nmod_poly_clear(_dr[i]);
	}
	nmod_poly_clear(tmp);
	nmod_poly_clear(d);
}

static int lin_verifier(nmod_poly_t y[WIDTH], nmod_poly_t _y[WIDTH],
		nmod_poly_t t, nmod_poly_t _t, nmod_poly_t u,
		commit_t com, commit_t x, commitkey_t *key,
		nmod_poly_t alpha, nmod_poly_t beta, int l,
		const uint8_t digest[SHA256HashSize]) {
	nmod_poly_t tmp, _d, v, _v, lin;
	int result = 1;

	nmod_poly_init(tmp, MODP);
	nmod_poly_init(lin, MODP);
	nmod_poly_init(_d, MODP);
	nmod_poly_init(v, MODP);
	nmod_poly_init(_v, MODP);
	nmod_poly_zero(v);
	nmod_poly_zero(_v);

	/* Sample challenge. */
	/* The proof carries the digest, not the first messages. */
	challenge_from_hash(_d, digest);

	/* Verifier checks norm. The responses are already ring elements, so there
	 * is nothing to reconstruct first. */
	for (int i = 0; i < WIDTH; i++) {
		/* Soundness checks: these must make verification fail rather than
		 * abort, and must not be compiled out by NDEBUG. */
		result &= commit_norm2_leq(y[i],
				(uint64_t) 4 * DEGREE * SIGMA_C * SIGMA_C);
		result &= commit_norm2_leq(_y[i],
				(uint64_t) 4 * DEGREE * SIGMA_C * SIGMA_C);
	}
	/* Verifier computes B1z and B1z'. */
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			commit_poly_mulmod(tmp, key->B1[i][j], y[j]);
			nmod_poly_add(v, v, tmp);
			commit_poly_mulmod(tmp, key->B1[i][j], _y[j]);
			nmod_poly_add(_v, _v, tmp);
		}
	}
	/* The first messages are not transmitted. Each is recovered from the
	 * equation that used to check it, t = B_1z - d c1, and the check happens
	 * once at the end when the digest is rebuilt over them. */
	commit_poly_mulmod(tmp, _d, com.c1);
	nmod_poly_sub(t, v, tmp);
	commit_poly_mulmod(tmp, _d, x.c1);
	nmod_poly_sub(_t, _v, tmp);

	/* The linear relation gets its own accumulator, because t now holds a
	 * recovered first message rather than a value about to be discarded. */
	commit_poly_mulmod(lin, alpha, com.c2);
	if (l == MSGS - 1 && (MSGS & 1)) {
		nmod_poly_sub(lin, lin, beta);
	} else {
		nmod_poly_add(lin, lin, beta);
	}
	nmod_poly_sub(lin, lin, x.c2);
	commit_poly_mulmod(lin, lin, _d);

	nmod_poly_zero(v);
	for (int i = 0; i < WIDTH; i++) {
		commit_poly_mulmod(tmp, key->b2[i], y[i]);
		commit_poly_mulmod(tmp, alpha, tmp);
		nmod_poly_add(v, v, tmp);
		commit_poly_mulmod(tmp, key->b2[i], _y[i]);
		nmod_poly_sub(v, v, tmp);
	}
	/* u = <b2 side> - d * <c2 side>, again a definition rather than a check. */
	nmod_poly_sub(u, v, lin);

	/* One comparison replaces the three first-message checks. */
	{
		uint8_t rebuilt[SHA256HashSize];
		nmod_poly_t ignored;

		nmod_poly_init(ignored, MODP);
		lin_hash(ignored, key, com, x, alpha, beta, u, t, _t, rebuilt);
		result &= (memcmp(rebuilt, digest, SHA256HashSize) == 0);
		nmod_poly_clear(ignored);
	}

	nmod_poly_clear(tmp);
	nmod_poly_clear(lin);
	nmod_poly_clear(_d);
	nmod_poly_clear(v);
	nmod_poly_clear(_v);
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
			hash_poly(&sha, _m[i]);
			/* No CRT splitting here: each component is a single polynomial. */
			hash_poly(&sha, com[i].c1);
			hash_poly(&sha, com[i].c2);
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
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_sub(sm[i], m[i], rho);
	}
	if (sc != NULL) {
		/* Ring elements are kept in coefficient representation, so rho is
		 * subtracted as it is, with no reduction into CRT components. */
		for (int i = 0; i < MSGS; i++) {
			nmod_poly_set(sc[i].c1, com[i].c1);
			nmod_poly_sub(sc[i].c2, com[i].c2, rho);
		}
	}
}

void shuffle_hash(nmod_poly_t beta, commit_t c[MSGS], commit_t d[MSGS],
		nmod_poly_t _m[MSGS], nmod_poly_t rho) {
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	prg_t prg;

	SHA256Reset(&sha);

	for (int i = 0; i < MSGS; i++) {
		SHA256Input(&sha, (const uint8_t *)_m[i]->coeffs,
				_m[i]->alloc * sizeof(uint64_t));
		SHA256Input(&sha, (const uint8_t *)c[i].c1->coeffs,
				c[i].c1->alloc * sizeof(uint64_t));
		SHA256Input(&sha, (const uint8_t *)c[i].c2->coeffs,
				c[i].c2->alloc * sizeof(uint64_t));
		SHA256Input(&sha, (const uint8_t *)d[i].c1->coeffs,
				d[i].c1->alloc * sizeof(uint64_t));
		SHA256Input(&sha, (const uint8_t *)d[i].c2->coeffs,
				d[i].c2->alloc * sizeof(uint64_t));
	}
	SHA256Input(&sha, (const uint8_t *)rho->coeffs,
			rho->alloc * sizeof(uint64_t));
	SHA256Result(&sha, hash);

	/* beta is public and must be reproducible by the verifier, so it is drawn
	 * from a stream keyed by the hash rather than from the global one. */
	prg_setseed(&prg, hash);
	commit_sample_rand(beta, &prg, DEGREE);
}

static void shuffle_prover(nmod_poly_t y[MSGS][WIDTH],
		nmod_poly_t _y[MSGS][WIDTH], nmod_poly_t t[MSGS],
		nmod_poly_t _t[MSGS], nmod_poly_t u[MSGS], commit_t d[MSGS],
		nmod_poly_t s[MSGS], commit_t com_in[MSGS], nmod_poly_t m_in[MSGS],
		nmod_poly_t _m_in[MSGS], nmod_poly_t r[MSGS][WIDTH],
		commitkey_t *key, uint8_t digest[MSGS][SHA256HashSize]) {
	nmod_poly_t beta, rho, t0, t1, theta[MSGS], _r[MSGS][WIDTH];
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
		for (int j = 0; j < WIDTH; j++) {
			nmod_poly_init(_r[i][j], MODP);
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
		commit_sample_short(_r[0][j]);
	}
	commit_doit(&d[0], t0, key, _r[0]);
	for (int i = 1; i < MSGS - 1; i++) {
		commit_sample_rand(theta[i], NULL, DEGREE);
		commit_poly_mulmod(t0, theta[i - 1], m[i]);
		commit_poly_mulmod(t1, theta[i], _m[i]);
		nmod_poly_add(t0, t0, t1);
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short(_r[i][j]);
		}
		commit_doit(&d[i], t0, key, _r[i]);
	}
	commit_poly_mulmod(t0, theta[MSGS - 2], m[MSGS - 1]);
	for (int j = 0; j < WIDTH; j++) {
		commit_sample_short(_r[MSGS - 1][j]);
	}
	commit_doit(&d[MSGS - 1], t0, key, _r[MSGS - 1]);

	shuffle_hash(beta, com, d, _m, rho);
	commit_poly_mulmod(s[0], theta[0], _m[0]);
	commit_poly_mulmod(t0, beta, m[0]);
	nmod_poly_sub(s[0], s[0], t0);
	nmod_poly_invmod(t0, _m[0], *commit_poly());
	commit_poly_mulmod(s[0], s[0], t0);
	for (int i = 1; i < MSGS - 1; i++) {
		commit_poly_mulmod(s[i], theta[i - 1], m[i]);
		commit_poly_mulmod(t0, theta[i], _m[i]);
		nmod_poly_add(s[i], s[i], t0);
		commit_poly_mulmod(t0, s[i - 1], m[i]);
		nmod_poly_sub(s[i], s[i], t0);
		nmod_poly_invmod(t0, _m[i], *commit_poly());
		commit_poly_mulmod(s[i], s[i], t0);
	}

	for (int l = 0; l < MSGS; l++) {
		if (l < MSGS - 1) {
			commit_poly_mulmod(t0, s[l], _m[l]);
		} else {
			commit_poly_mulmod(t0, beta, _m[l]);
		}

		if (l == 0) {
			lin_prover(digest[l], y[l], _y[l], t[l], _t[l], u[l], com[l], d[l], key, beta,
					t0, r[l], _r[l], l);
		} else {
			lin_prover(digest[l], y[l], _y[l], t[l], _t[l], u[l], com[l], d[l], key,
					s[l - 1], t0, r[l], _r[l], l);
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
		for (int j = 0; j < WIDTH; j++) {
			nmod_poly_clear(_r[i][j]);
		}
	}
}

static int shuffle_verifier(nmod_poly_t y[MSGS][WIDTH],
		nmod_poly_t _y[MSGS][WIDTH], nmod_poly_t t[MSGS],
		nmod_poly_t _t[MSGS], nmod_poly_t u[MSGS], commit_t d[MSGS],
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
 * the shuffled list _m[] are the statement, and the key is a public parameter.
 * Neither are the first messages t, _t and u: the verifier recovers those from
 * the equations, and the proof carries the 32-byte digest in their place.
 */
typedef struct _proof_t {
	nmod_poly_t (*y)[WIDTH];
	nmod_poly_t (*_y)[WIDTH];
	commit_t *d;
	nmod_poly_t *s;
	uint8_t (*digest)[SHA256HashSize];
} proof_t;

static void walk_uniform(bitwriter_t *w, bitreader_t *r, nmod_poly_t a) {
	if (w != NULL) {
		serial_put_uniform(w, a);
	} else {
		serial_get_uniform(r, a);
	}
}

static void walk_gauss(bitwriter_t *w, bitreader_t *r, nmod_poly_t a,
		ulong sigma) {
	if (w != NULL) {
		serial_put_gauss(w, a, sigma);
	} else {
		serial_get_gauss(r, a, sigma);
	}
}

static void proof_walk(proof_t *pf, bitwriter_t *w, bitreader_t *r) {
	for (int l = 0; l < MSGS; l++) {
		for (int i = 0; i < WIDTH; i++) {
			walk_gauss(w, r, pf->y[l][i], SIGMA_C);
			walk_gauss(w, r, pf->_y[l][i], SIGMA_C);
		}
		walk_uniform(w, r, pf->d[l].c1);
		walk_uniform(w, r, pf->d[l].c2);
		walk_uniform(w, r, pf->s[l]);
		for (int i = 0; i < SHA256HashSize; i++) {
			if (w != NULL) {
				serial_put_byte(w, pf->digest[l][i]);
			} else {
				pf->digest[l][i] = serial_get_byte(r);
			}
		}
	}
}

/* Set by the size test: round-trip the proof through its serialisation before
 * verifying, and record how many bytes it took. */
static int measure_proof = 0;
static size_t proof_bytes = 0;

static int run(commit_t com[MSGS], nmod_poly_t m[MSGS], nmod_poly_t _m[MSGS],
		nmod_poly_t r[MSGS][WIDTH], commitkey_t *key) {
	int result = 1;
	commit_t d[MSGS];
	nmod_poly_t s[MSGS], u[MSGS];
	nmod_poly_t y[MSGS][WIDTH], _y[MSGS][WIDTH], t[MSGS], _t[MSGS];

	for (int i = 0; i < MSGS; i++) {
		commit_init(&d[i]);
		nmod_poly_init(s[i], MODP);
		nmod_poly_init(t[i], MODP);
		nmod_poly_init(_t[i], MODP);
		nmod_poly_init(u[i], MODP);
		for (int j = 0; j < WIDTH; j++) {
			nmod_poly_init(y[i][j], MODP);
			nmod_poly_init(_y[i][j], MODP);
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
		serial_reader_init(&br, buf, sizeof(buf));
		proof_walk(&pf, NULL, &br);
	}

	result = shuffle_verifier(y, _y, t, _t, u, d, s, com, _m, key, digest);

	for (int i = 0; i < MSGS; i++) {
		commit_free(&d[i]);
		nmod_poly_clear(s[i]);
		nmod_poly_clear(t[i]);
		nmod_poly_clear(_t[i]);
		nmod_poly_clear(u[i]);
		for (int j = 0; j < WIDTH; j++) {
			nmod_poly_clear(y[i][j]);
			nmod_poly_clear(_y[i][j]);
		}
	}

	return result;
}

static void test(void) {
	commitkey_t key;
	commit_t com[MSGS], com0[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS], m0[MSGS], _m0[MSGS], r[MSGS][WIDTH];

	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
		nmod_poly_init(m0[i], MODP);
		nmod_poly_init(_m0[i], MODP);
		commit_init(&com0[i]);
		for (int j = 0; j < WIDTH; j++) {
			nmod_poly_init(r[i][j], MODP);
		}
	}

	/* Generate commitment key-> */
	commit_setup();
	commit_keyinit(&key);
	commit_keygen(&key);

	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short(r[i][j]);
		}
		commit_sample_short(m[i]);
		commit_init(&com[i]);
		commit_doit(&com[i], m[i], &key, r[i]);
	}

	/* Prover shuffles messages (only a circular shift for simplicity). */
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_set(_m[i], m[(i + 1) % MSGS]);
	}

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
			double fixed = 3.0 * DEGREE * serial_bits_uniform()
					+ SHA256HashSize * 8;
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
			nmod_poly_set(com0[i].c1, com[i].c1);
			nmod_poly_set(com0[i].c2, com[i].c2);
		}
		TEST_ASSERT(run(com, m, _m, r, &key) == 1, end);
		for (int i = 0; i < MSGS; i++) {
			TEST_ASSERT(nmod_poly_equal(m0[i], m[i]) == 1, end);
			TEST_ASSERT(nmod_poly_equal(_m0[i], _m[i]) == 1, end);
			TEST_ASSERT(nmod_poly_equal(com0[i].c1, com[i].c1), end);
			TEST_ASSERT(nmod_poly_equal(com0[i].c2, com[i].c2), end);
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
			nmod_poly_clear(r[i][j]);
		}
	}
	commit_keyfree(&key);
}

static void bench(void) {
	commitkey_t key;
	commit_t com[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS];
	nmod_poly_t alpha, beta, s[MSGS - 1];
	nmod_poly_t r[MSGS][WIDTH], y[WIDTH], _y[WIDTH];
	nmod_poly_t t, _t, u, v, _v;

	nmod_poly_init(alpha, MODP);
	nmod_poly_init(beta, MODP);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
		if (i != MSGS - 1)
			nmod_poly_init(s[i], MODP);
		for (int j = 0; j < WIDTH; j++) {
			nmod_poly_init(r[i][j], MODP);
		}
	}

	/* Generate commitment key-> */
	commit_setup();
	commit_keyinit(&key);
	commit_keygen(&key);

	for (int i = 0; i < MSGS; i++) {
		for (int j = 0; j < WIDTH; j++) {
			commit_sample_short(r[i][j]);
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

	nmod_poly_init(t, MODP);
	nmod_poly_init(_t, MODP);
	nmod_poly_init(u, MODP);
	nmod_poly_init(v, MODP);
	nmod_poly_init(_v, MODP);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_init(y[i], MODP);
		nmod_poly_init(_y[i], MODP);
	}

	/* The randomness must stay the one the commitments were formed with:
	 * resampling it here would benchmark a proof whose witness does not match
	 * its statement, and a verifier that always rejects. */
	commit_sample_rand(beta, NULL, DEGREE);
	commit_sample_rand(alpha, NULL, DEGREE);
	uint8_t dg[SHA256HashSize];

	BENCH_BEGIN("linear proof") {
		BENCH_ADD(lin_prover(dg, y, _y, t, _t, u, com[0], com[1], &key, alpha, beta,
						r[0], r[0], 0));
	} BENCH_END;

	BENCH_BEGIN("linear verifier") {
		BENCH_ADD(lin_verifier(y, _y, t, _t, u, com[0], com[1], &key, alpha,
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
			nmod_poly_clear(r[i][j]);
		}
	}
	nmod_poly_clear(t);
	nmod_poly_clear(_t);
	nmod_poly_clear(u);
	nmod_poly_clear(v);
	nmod_poly_clear(_v);
	for (int i = 0; i < WIDTH; i++) {
		nmod_poly_clear(y[i]);
		nmod_poly_clear(_y[i]);
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
