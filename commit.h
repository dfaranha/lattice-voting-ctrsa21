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

#include <flint/flint.h>
#include <flint/nmod_poly.h>

#include "param.h"

/*============================================================================*/
/* Constant definitions                                                       */
/*============================================================================*/

/* Parameter v in the commitment scheme (laximum l1-norm of challs). */
#define NONZERO 36
/* The \infty-norm bound of certain elements. */
#define BETA 	1
/* Width k of the comming matrix. */
#define WIDTH 	3
/* Height of the commitment matrix. */
#define HEIGHT 	1

/*============================================================================*/
/* Type definitions                                                           */
/*============================================================================*/

/* Type that represents a commitment key pair. Ring elements are kept in
 * coefficient representation: this build does not exploit the CRT splitting of
 * the ring for arithmetic. */
typedef struct _key_t {
	nmod_poly_t B1[HEIGHT][WIDTH];
	nmod_poly_t b2[WIDTH];
} commitkey_t;

/* Type that represents a commitment. */
typedef struct _com_t {
	nmod_poly_t c1, c2;
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

/* Multiply two polynomials in the cyclotomic ring Rp.
 *
 * @param[out] c		- the resulting polynomial.
 * @param[in] a			- the first polynomial.
 * @param[in] b			- the second polynomial.
 */
void commit_poly_mulmod(nmod_poly_t c, const nmod_poly_t a,
		const nmod_poly_t b);

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
 * Sample a random polynomial.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_rand(nmod_poly_t r, flint_rand_t rand, int degree);

/**
 * Sample a random challenge.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_chall(nmod_poly_t f);

/**
 * Sample a random polynomial following a Gaussian distribution.
 *
 * @param[out] r		- the polynomial to sample.
 */
void commit_sample_gauss(nmod_poly_t r);

/**
 * Initialize a commitment. Must be called before commit_doit, and released
 * with commit_free. Separating this from commit_doit is what allows the same
 * commitment to be recomputed in a loop without leaking.
 *
 * @param[out] com 		- the commitment to initialize.
 */
void commit_init(commit_t *com);

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
		nmod_poly_t r[WIDTH]);

/**
 * Open a commitment on a certain message.
 *
 * @param[in] com 		- the commitment to open.
 * @param[in] m 		- the associated message.
 * @param[in] r			- the opening randomness.
 * @param[in] f			- the opening challenge.
 */
int commit_open(commit_t *com, nmod_poly_t m, commitkey_t *key,
		nmod_poly_t r[WIDTH], nmod_poly_t f);

/**
 * Free a commitment.
 *
 * @param[out] com 		- the commitment to free.
 */
void commit_free(commit_t *com);
