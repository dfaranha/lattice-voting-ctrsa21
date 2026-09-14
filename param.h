/**
 * @defgroup evoting Lattice-based electronic voting.
 */
/**
 * @file
 *
 * Common parameters to all lattice-based schemes.
 *
 * @ingroup evoting
 */

#ifndef PARAM_H
#define PARAM_H

#include <sys/random.h>
#include <stdint.h>
#include <flint/flint.h>
#include <flint/fmpz_mod_poly.h>
#include <string.h>

/* FLINT 3 renamed these, deprecating the old spellings. Map the new names onto
 * the old ones when building against an earlier release, so that the sources
 * can use the current API either way. */
#if !defined(__FLINT_RELEASE) || __FLINT_RELEASE < 30000
#define flint_rand_init 	flint_randinit
#define flint_rand_clear 	flint_randclear
#define flint_rand_set_seed 	flint_randseed
#endif

/*============================================================================*/
/* Constant definitions                                                       */
/*============================================================================*/

/* Modulus p defining the cyclotomic ring for the commitment scheme. */
#define MODP 	3906450253
/* Degree of the polynomial defining the cyclotomic ring for the commitment scheme. */
#ifndef DEGREE
#define DEGREE 	1024
#endif
/* Degree of each polynomial used to define the CRT representation. */
#define DEGCRT 	(DEGREE >> 1)
/* Number of irreducible factors of (x^DEGREE + 1), and therefore the number of
 * components of a CRT representation. This is a property of the ring, and is
 * unrelated to the module dimensions DIM, WIDTH and HEIGHT: it only happens to
 * share the value 2 with DIM, which made several loops index a CRT component
 * with a module bound. */
#define NCRT 	2
/* Standard deviation for discrete Gaussians. */
#define SIGMA_C 54000
/* Standard deviation for the Gaussian masking the committed permutation
 * elements sigma_i of Lemma 5. These are monomials, so ||d * sigma_i|| equals
 * ||d|| <= sqrt(2 * NONZERO) ~ 8.5 exactly, and SIGMA_S = 256 leaves a ratio of
 * about 30 for the rejection sampling. It must also stay small enough that the
 * extracted ||(d - d') * sigma_i|| <= 2 * (2 * sqrt(DEGREE) * SIGMA_S) = 32768
 * remains below sqrt(MODP) ~ 62501, the bound under which Lemma 1 guarantees
 * invertibility. */
#define SIGMA_S 256
// Standard deviation of discrete Gaussian
#define SIGMA_E	54000

#endif /* PARAM_H */
