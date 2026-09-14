/**
 * @defgroup commit Lattice-based commitment.
 */
/**
 * @file
 *
 * Interface of the lattice-based commitment scheme.
 *
 * @ingroup commit
 */

#ifndef COMMIT_H
#define COMMIT_H

#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>

#include "param.h"

/*============================================================================*/
/* Constant definitions                                                       */
/*============================================================================*/

/* Parameter v in the commitment scheme (laximum l1-norm of challs). */
#define NONZERO 36
/* The \infty-norm bound of certain elements. */
#define BETA 	1
/* Width k of the commitment matrix. The hiding property is MLWE of rank
 * k - n - 1, which at k = 3 is rank 1 and only 49.6 bits at this modulus.
 * k = 4 gives rank 2. */
#define WIDTH 	4
/* Height of the commitment matrix. */
#define HEIGHT 	1

/*============================================================================*/
/* Type definitions                                                           */
/*============================================================================*/

/* Type that represents a polynomial in CRT representation. */
typedef nmod_poly_t pcrt_poly_t[NCRT];

/* Type that represents a commitment key pair. */
typedef struct _key_t {
	pcrt_poly_t B1[HEIGHT][WIDTH];
	pcrt_poly_t b2[WIDTH];
} commitkey_t;

/* Type that represents a commitment in CRT representation. */
typedef struct _com_t {
	pcrt_poly_t c1, c2;
} commit_t;

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Initialize the commitment module.
 */
void commit_setup();

/**
 * Finalize the commitment module.
 */
void commit_finish();

/**
 * Return the polynomial (x^N + 1) defining the cyclotomic ring Rp.
 *
 * @return the polynomial defining the cyclotomic ring Rp.
 */
nmod_poly_t *commit_poly();

/**
 * Return the i-th polynomial defining the CRT representation.
 * The returned polynomial is such that it factors the polynomial defining Rp.
 *
 * @return the i-th irreducible polynomial defining the CRT representation.
 */
nmod_poly_t *commit_irred(int i);

/* Multiply two polynomials modulo the i-th CRT factor.
 *
 * @param[out] c		- the resulting polynomial.
 * @param[in] a			- the first polynomial.
 * @param[in] b			- the second polynomial.
 * @param[in] i			- the index of the CRT factor.
 */
void pcrt_poly_mulmod(nmod_poly_t c, const nmod_poly_t a, const nmod_poly_t b,
		int i);

/* Multiply two polynomials in the cyclotomic ring Rp.
 *
 * @param[out] c		- the resulting polynomial.
 * @param[in] a			- the first polynomial.
 * @param[in] b			- the second polynomial.
 */
void commit_poly_mulmod(nmod_poly_t c, const nmod_poly_t a,
		const nmod_poly_t b);

/* Reduce a polynomial into the i-th CRT component.
 *
 * The output may alias the input.
 *
 * @param[out] c		- the reduced polynomial.
 * @param[in] a			- the polynomial to reduce.
 * @param[in] i			- the index of the CRT factor.
 */
void pcrt_poly_reduce(nmod_poly_t c, const nmod_poly_t a, int i);

/* Recover polynomial from CRT representation.
 *
 * @param[in] c 		- the resulting polynomial.
 * @param[in] a 		- the polynomial in CRT representation.
 */
void pcrt_poly_rec(nmod_poly_t c, pcrt_poly_t a);

/* Compute the squared l2-norm of a polynomial.
 *
 * @param[in] r			- the polynomial to compute the norm.
 * @return The squared l2-norm.
 */
uint64_t commit_norm2_sqr(nmod_poly_t r);

/* Compute the l\infty-norm of a polynomial.
 *
 * @param[in] r			- the polynomial to compute the norm.
 * @return The l\infty-norm.
 */
uint64_t commit_norm_inf(nmod_poly_t r);

/* Test whether the squared l2-norm of a polynomial is at most a bound.
 *
 * Unlike commit_norm2_sqr this cannot overflow for any modulus, so it is safe
 * on input chosen
 * by a malicious party.
 *
 * @param[in] r			- the polynomial to test.
 * @param[in] bound		- the bound on the squared l2-norm.
 * @return 1 if the squared l2-norm is at most the bound, 0 otherwise.
 */
int commit_norm2_leq(nmod_poly_t r, uint64_t bound);

/**
 * Initialize a key pair for the commitment scheme. Must be called before
 * commit_keygen, and released with commit_keyfree.
 *
 * @param[out] key 		- the key pair to initialize.
 */
void commit_keyinit(commitkey_t *key);

/**
 * Generate a key pair for the commitment scheme using a PRNG.
 *
 * @param[out] key 		- the generated key pair.
 * @param[in] rand		- the random number generator.
 */
void commit_keygen(commitkey_t *key, flint_rand_t rand);

/**
 * Free a key pair for the commitment scheme.
 *
 * @param[out] key 		- the generated key pair.
 * @param[in] rand		- the random number generator.
 */
void commit_keyfree(commitkey_t *key);

/**
 * Sample a short polynomial.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_short(nmod_poly_t r);

/**
 * Sample a short polynomial in CRT representation.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_short_crt(pcrt_poly_t r);

/**
 * Sample a random polynomial.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_rand(nmod_poly_t r, flint_rand_t rand, int degree);

/**
 * Sample a random polynomial in CRT representation.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_rand_crt(pcrt_poly_t r, flint_rand_t rand);

/**
 * Sample a random challenge.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_chall(nmod_poly_t f);

/**
 * Sample a random challenge in CRT representation.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_chall_crt(pcrt_poly_t f);

/**
 * Sample a random polynomial following a Gaussian distribution.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_gauss(nmod_poly_t r);

/**
 * Sample a random polynomial following a Gaussian distribution in CRT
 * representation.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_gauss_crt(pcrt_poly_t r);

/**
 * Sample a masking polynomial at the batch width SIGMA_B, for a proof whose
 * rejection test spans every message at once.
 */
void commit_sample_gauss_batch_crt(pcrt_poly_t r);

/**
 * Initialize a commitment. Must be called before commit_doit, and released
 * with commit_free. Separating this from commit_doit is what allows the same
 * commitment to be recomputed in a loop without leaking.
 *
 * @param[out] com 		- the commitment to initialize.
 */
void commit_init(commit_t *com);

/**
 * Sample a random polynomial following a narrow Gaussian distribution of
 * standard deviation SIGMA_S, in CRT representation. Used to mask committed
 * messages that must be proven short.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_gauss_small_crt(pcrt_poly_t r);

/**
 * Commit to a message and randomness using a key pair.
 *
 * The commitment must already be initialized with commit_init.
 *
 * @param[out] com 		- the resulting commitment.
 * @param[in] m 		- the message to commit.
 * @param[in] r 		- the commitment randomness.
 */
void commit_doit(commit_t *com, nmod_poly_t m, commitkey_t *key,
		pcrt_poly_t r[WIDTH]);

/**
 * Open a commitment on a certain message.
 *
 * @param[in] com 		- the commitment to open.
 * @param[in] m 		- the associated message.
 * @param[in] r			- the opening randomness.
 * @param[in] f			- the opening challenge.
 */
int commit_open(commit_t *com, nmod_poly_t m, commitkey_t *key,
		pcrt_poly_t r[WIDTH], pcrt_poly_t f);

/**
 * Free a commitment.
 *
 * @param[out] com 		- the commitment to free.
 */
void commit_free(commit_t *com);

/**
 * Sample a uniformly random double in [0, 1) with 53 bits of precision.
 *
 * Rejection sampling must not reuse the Fiat-Shamir stream, so this draws from
 * the operating system rather than from fastrandombytes.
 *
 * @return the sampled value.
 */
double commit_uniform_double(void);

/**
 * Decide whether a masked response must be rejected, following Fiat-Shamir
 * with aborts. Shared by the shuffle proof and the LNP machinery so that this
 * security-critical step exists only once.
 *
 * @param[in] z				- the masked responses, in CRT representation.
 * @param[in] v				- the terms they mask, in CRT representation.
 * @param[in] s2			- the square of the standard deviation.
 * @param[in] width			- the number of responses.
 * @return 1 if the transcript must be rejected, 0 if it may be emitted.
 */
int commit_rej_sampling(nmod_poly_t z[][2], nmod_poly_t v[][2], uint64_t s2,
		int width);

/**
 * Accumulate the two quantities the rejection test needs, so that one decision
 * can span responses held in several separate arrays. A proof that shares a
 * challenge across many statements has to test their concatenation, not each
 * piece: testing them separately would multiply the abort probabilities.
 *
 * @param[in,out] dot		- running inner product of the responses and the
 *							  terms they mask.
 * @param[in,out] norm		- running squared norm of the masked terms.
 */
void commit_rej_accumulate(int64_t *dot, int64_t *norm, nmod_poly_t z[][2],
		nmod_poly_t v[][2], int width);

/**
 * Decide from accumulated quantities.
 *
 * @return 1 if the transcript must be rejected, 0 if it may be emitted.
 */
int commit_rej_decide(int64_t dot, int64_t norm, uint64_t s2);

#endif /* COMMIT_H */
