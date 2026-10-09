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
#include "fastrandombytes.h"
#include "sha.h"

/*============================================================================*/
/* Constant definitions                                                       */
/*============================================================================*/

/* Parameter nu in the commitment scheme: the number of non-zero coefficients,
 * each +-1, of a challenge, and so its l1-norm. A difference of two
 * challenges, which is what an opening carries, has l1-norm at most
 * 2 * NONZERO. The challenge space has 2^NONZERO * C(DEGREE, NONZERO)
 * elements, which at DEGREE 2048 is 2^163.4 even at 18, comfortably past the
 * 2^128 the knowledge error needs. Lowering it pays twice: SIGMA_C is linear
 * in NONZERO, so the masks narrow, and the extractable opening
 * 16 sigma sqrt(nu N) falls with both, taking MSIS binding from 309 to 352
 * bits. */
#define NONZERO 18
/* The \infty-norm bound of certain elements. */
#define BETA 	1
/* Width k of the commitment matrix. The hiding property is MLWE of rank
 * k - n - 1, so k = 3 gives rank 1 and k = 4 rank 2. Rank 1 was not an option
 * at DEGREE 1024, where it is about 52 bits, but at 2048 it is 140 and rank 2
 * would be 343, far past what is needed. So k = 3 here. */
#define WIDTH 	3
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
	/* A hash of B1 and b2, set by commit_keygen. Every Fiat-Shamir challenge
	 * binds the key, and absorbing this instead of the matrices themselves
	 * saves rehashing a constant on every attempt of every proof. */
	uint8_t digest[SHA256HashSize];
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

/* Multiply two polynomials in the cyclotomic ring Rp, the first of them sparse.
 *
 * The result is the same as commit_poly_mulmod's, at a cost proportional to the
 * number of non-zero coefficients of d, which makes it the product to use for
 * a challenge. The output may alias either input.
 *
 * @param[out] c		- the resulting polynomial.
 * @param[in] d			- the sparse polynomial.
 * @param[in] a			- the other polynomial.
 */
void commit_poly_mul_sparse(nmod_poly_t c, const nmod_poly_t d,
		const nmod_poly_t a);

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
 * Generate a key pair for the commitment scheme.
 *
 * @param[out] key 		- the generated key pair.
 */
void commit_keygen(commitkey_t *key);

/**
 * Absorb a polynomial into a hash computation in a canonical way.
 *
 * The coefficients are absorbed one by one up to DEGREE, so that the digest
 * depends only on the value of the polynomial and not on the capacity FLINT
 * happens to have allocated. Each takes the fewest bytes that hold MODP: four
 * when MODP < 2^32, eight otherwise.
 *
 * @param[in,out] sha		- the hash context.
 * @param[in] p				- the polynomial to absorb.
 */
void commit_hash_poly(SHA256Context *sha, nmod_poly_t p);

/**
 * Absorb a CRT component into a hash computation in a canonical way.
 *
 * A component is reduced modulo a factor of degree DEGCRT, so only its first
 * DEGCRT coefficients can be non-zero, and only those are absorbed. Absorbing
 * DEGREE of them, as hash_poly used to, hashed half a component of zeros.
 * A component that is not reduced is a bug, and aborts rather than hashing
 * a truncation that two different values could share.
 *
 * @param[in,out] sha		- the hash context.
 * @param[in] p				- the component to absorb.
 */
void commit_hash_crt(SHA256Context *sha, nmod_poly_t p);

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
 * Sample a uniformly random polynomial.
 *
 * @param[out] r		- the polynomial to sample.
 * @param[in] prg		- the stream to draw from, or NULL for the global one,
 * 						  which is seeded by the operating system and is what
 * 						  anything secret must use.
 * @param[in] degree	- the number of coefficients to sample.
 */
void commit_sample_rand(nmod_poly_t r, prg_t *prg, int degree);

/**
 * Sample a uniformly random polynomial in CRT representation.
 *
 * @param[out] r		- the polynomial to sample.
 * @param[in] prg		- the stream to draw from, or NULL for the global one.
 */
void commit_sample_rand_crt(pcrt_poly_t r, prg_t *prg);

/**
 * Sample a random opening challenge: the difference of two challenges with
 * NONZERO coefficients of +-1 each, which is what an extracted opening
 * carries.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_chall(nmod_poly_t f);

/**
 * Sample a random opening challenge in CRT representation.
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
 *
 * @param[out] r		- the polynomial to sample, in coefficient form.
 */
void commit_sample_gauss_batch(nmod_poly_t r);

/**
 * Sample a masking polynomial at the batch width SIGMA_B, in CRT
 * representation.
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
 * @param[in] z				- the masked responses, in coefficient form.
 * @param[in] v				- the terms they mask, in coefficient form.
 * @param[in] s2			- the square of the standard deviation.
 * @param[in] width			- the number of responses.
 * @return 1 if the transcript must be rejected, 0 if it may be emitted.
 */
int commit_rej_sampling(nmod_poly_t z[], nmod_poly_t v[], uint64_t s2,
		int width);

/**
 * Accumulate the two quantities the rejection test needs, so that one decision
 * can span responses held in several separate arrays. A proof that shares a
 * challenge across many statements has to test their concatenation, not each
 * piece: testing them separately would multiply the abort probabilities.
 *
 * The responses are taken in coefficient form, which is what the test reads:
 * reconstructing them from CRT form costs a product at full degree each, and
 * a prover that keeps its masks and products in coefficient form need not pay
 * it on every attempt.
 *
 * @param[in,out] dot		- running inner product of the responses and the
 *							  terms they mask.
 * @param[in,out] norm		- running squared norm of the masked terms.
 * @param[in] z				- the masked responses, in coefficient form.
 * @param[in] v				- the terms they mask, in coefficient form.
 * @param[in] width			- the number of responses.
 */
void commit_rej_accumulate(int64_t *dot, int64_t *norm, nmod_poly_t z[],
		nmod_poly_t v[], int width);

/**
 * Answer a challenge: z = y + d * v for each of width masks y and witnesses v,
 * accumulating what the rejection test needs and reducing each response into
 * the CRT form a proof carries. Everything is in coefficient form except the
 * output, so that the test reads the responses without reconstructing them,
 * and d being a challenge, each product is a sparse rotation.
 *
 * @param[out] z			- the responses, in CRT representation.
 * @param[in,out] y			- the masks; on return, the responses, both in
 *							  coefficient form.
 * @param[in] v				- the witness, in coefficient form.
 * @param[in] d				- the challenge, in coefficient form.
 * @param[in] width			- the number of responses.
 * @param[in,out] dot		- as in commit_rej_accumulate.
 * @param[in,out] norm		- as in commit_rej_accumulate.
 */
void commit_rej_respond(nmod_poly_t z[][2], nmod_poly_t y[], nmod_poly_t v[],
		nmod_poly_t d, int width, int64_t *dot, int64_t *norm);

/**
 * Decide from accumulated quantities.
 *
 * @return 1 if the transcript must be rejected, 0 if it may be emitted.
 */
int commit_rej_decide(int64_t dot, int64_t norm, uint64_t s2);

#endif /* COMMIT_H */
