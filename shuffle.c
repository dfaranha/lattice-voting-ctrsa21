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
	uint32_t buf;

	fastrandombytes_setseed((uint8_t *) hash);
	/* The challenge is the difference of two polynomials with NONZERO ones.
	 * The 2 here belongs to that construction, not to any CRT splitting. */
	nmod_poly_t c[2];
	for (int i = 0; i < 2; i++) {
		nmod_poly_init(c[i], MODP);
		nmod_poly_fit_length(c[i], DEGREE);
		for (int j = 0; j < NONZERO; j++) {
			fastrandombytes((unsigned char *)&buf, sizeof(buf));
			buf = buf % DEGREE;
			while (nmod_poly_get_coeff_ui(c[i], buf) != 0) {
				fastrandombytes((unsigned char *)&buf, sizeof(buf));
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
	// Compute sigma^2 = (11 * v * beta * sqrt(k * N))^2.
	uint64_t sigma_sqr = 11 * NONZERO * BETA;
	sigma_sqr *= sigma_sqr * DEGREE * WIDTH;

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

void shuffle_hash(nmod_poly_t beta, commit_t c[MSGS], commit_t d[MSGS],
		nmod_poly_t _m[MSGS], nmod_poly_t rho) {
	flint_rand_t rand;
	SHA256Context sha;
	uint8_t hash[SHA256HashSize];
	uint64_t seed0, seed1, seed2, seed3;

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

	flint_rand_init(rand);
	memcpy(&seed0, hash, sizeof(uint64_t));
	memcpy(&seed1, hash + sizeof(uint64_t), sizeof(uint64_t));
	memcpy(&seed2, hash + 2 * sizeof(uint64_t), sizeof(uint64_t));
	memcpy(&seed3, hash + 3 * sizeof(uint64_t), sizeof(uint64_t));
	seed0 ^= seed2;
	seed1 ^= seed3;
	flint_rand_set_seed(rand, seed0, seed1);
	commit_sample_rand(beta, rand, DEGREE);
	flint_rand_clear(rand);
}

static void shuffle_prover(nmod_poly_t y[MSGS][WIDTH],
		nmod_poly_t _y[MSGS][WIDTH], nmod_poly_t t[MSGS],
		nmod_poly_t _t[MSGS], nmod_poly_t u[MSGS], commit_t d[MSGS],
		nmod_poly_t s[MSGS], commit_t com[MSGS], nmod_poly_t m[MSGS],
		nmod_poly_t _m[MSGS], nmod_poly_t r[MSGS][WIDTH], nmod_poly_t rho,
		commitkey_t *key, uint8_t digest[MSGS][SHA256HashSize],
		flint_rand_t rng) {
	nmod_poly_t beta, t0, t1, theta[MSGS], _r[MSGS][WIDTH];

	nmod_poly_init(t0, MODP);
	nmod_poly_init(t1, MODP);
	nmod_poly_init(beta, MODP);
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(theta[i], MODP);
		for (int j = 0; j < WIDTH; j++) {
			nmod_poly_init(_r[i][j], MODP);
		}
	}

	/* Prover shifts the messages by rho. */
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_sub(m[i], m[i], rho);
		nmod_poly_sub(_m[i], _m[i], rho);
	}

	/* Prover samples theta_i and computes commitments D_i. */
	commit_sample_rand(theta[0], rng, DEGREE);
	commit_poly_mulmod(t0, theta[0], _m[0]);
	for (int j = 0; j < WIDTH; j++) {
		commit_sample_short(_r[0][j]);
	}
	commit_doit(&d[0], t0, key, _r[0]);
	for (int i = 1; i < MSGS - 1; i++) {
		commit_sample_rand(theta[i], rng, DEGREE);
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
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_clear(theta[i]);
		for (int j = 0; j < WIDTH; j++) {
			nmod_poly_clear(_r[i][j]);
		}
	}
}

static int shuffle_verifier(nmod_poly_t y[MSGS][WIDTH],
		nmod_poly_t _y[MSGS][WIDTH], nmod_poly_t t[MSGS],
		nmod_poly_t _t[MSGS], nmod_poly_t u[MSGS], commit_t d[MSGS],
		nmod_poly_t s[MSGS], commit_t com[MSGS], nmod_poly_t _m[MSGS],
		nmod_poly_t rho, commitkey_t *key,
		uint8_t digest[MSGS][SHA256HashSize]) {
	int result = 1;
	nmod_poly_t beta, t0;

	nmod_poly_init(t0, MODP);
	nmod_poly_init(beta, MODP);

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
		nmod_poly_t r[MSGS][WIDTH], commitkey_t *key, flint_rand_t rng) {
	int flag, result = 1;
	commit_t d[MSGS];
	nmod_poly_t t0, t1, rho, s[MSGS], u[MSGS];
	nmod_poly_t y[MSGS][WIDTH], _y[MSGS][WIDTH], t[MSGS], _t[MSGS];

	nmod_poly_init(t0, MODP);
	nmod_poly_init(t1, MODP);
	nmod_poly_init(rho, MODP);
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

	/* Verifier samples \rho that is different from the messages, and \beta. */
	do {
		flag = 1;
		commit_sample_rand(rho, rng, DEGREE);
		for (int i = 0; i < MSGS; i++) {
			if (nmod_poly_equal(rho, _m[i]) == 1) {
				flag = 0;
			}
		}
	} while (flag == 0);

	/* Verifier shifts the commitments by rho. */
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_sub(com[i].c2, com[i].c2, rho);
	}

	static uint8_t digest[MSGS][SHA256HashSize];

	shuffle_prover(y, _y, t, _t, u, d, s, com, m, _m, r, rho, key, digest, rng);

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

	result = shuffle_verifier(y, _y, t, _t, u, d, s, com, _m, rho, key,
			digest);

	nmod_poly_clear(t0);
	nmod_poly_clear(t1);
	nmod_poly_clear(rho);
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

static void test(flint_rand_t rand) {
	commitkey_t key;
	commit_t com[MSGS];
	nmod_poly_t m[MSGS], _m[MSGS], r[MSGS][WIDTH];

	for (int i = 0; i < MSGS; i++) {
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
		for (int j = 0; j < WIDTH; j++) {
			nmod_poly_init(r[i][j], MODP);
		}
	}

	/* Generate commitment key-> */
	commit_setup();
	commit_keyinit(&key);
	commit_keygen(&key, rand);

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
		TEST_ASSERT(run(com, m, _m, r, &key, rand) == 1, end);
		measure_proof = 0;
		printf("\n    %zu bytes, %.1f KB per message\n", proof_bytes,
				proof_bytes / 1024.0 / MSGS);
		{
			size_t bu = serial_bits_uniform();
			size_t bg = serial_bits_gauss(SIGMA_C);
			size_t per = 3 * DEGREE * bu + 2 * WIDTH * DEGREE * bg
					+ SHA256HashSize * 8;

			TEST_ASSERT(proof_bytes == (MSGS * per + 7) / 8, end);
		}
	} TEST_END;

	TEST_ONCE("shuffle proof is consistent") {
		TEST_ASSERT(run(com, m, _m, r, &key, rand) == 1, end);
	} TEST_END;

  end:
	commit_finish();

	for (int i = 0; i < MSGS; i++) {
		commit_free(&com[i]);
		nmod_poly_clear(m[i]);
		nmod_poly_clear(_m[i]);
		for (int j = 0; j < WIDTH; j++) {
			nmod_poly_clear(r[i][j]);
		}
	}
	commit_keyfree(&key);
}

static void bench(flint_rand_t rand) {
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
	commit_keygen(&key, rand);

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
		BENCH_ADD(run(com, m, _m, r, &key, rand));
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
	commit_sample_rand(beta, rand, DEGREE);
	commit_sample_rand(alpha, rand, DEGREE);
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
