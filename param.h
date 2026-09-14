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

 #include <sys/random.h>
#include <stdint.h>
#include <flint/flint.h>
#include <flint/fmpz_mod_poly.h>
#include <string.h>

/*============================================================================*/
/* Constant definitions                                                       */
/*============================================================================*/

/* Modulus p defining the cyclotomic ring for the commitment scheme. */
#define MODP 	3906450253
/* Degree of the polynomial defining the cyclotomic ring for the commitment scheme. */
#define DEGREE 	1024
/* Degree of each polynomial used to define the CRT representation. */
#define DEGCRT 	(DEGREE >> 1)
/* Standard deviation for discrete Gaussians. */
#define SIGMA_C 54000
// Standard deviation of discrete Gaussian
#define SIGMA_E	54000
/* Standard deviation for the Gaussian masking the committed permutation
 * elements sigma_i of Lemma 5. These are monomials, so ||d * sigma_i|| equals
 * ||d|| <= sqrt(2 * NONZERO) ~ 8.5 exactly, and SIGMA_S = 256 leaves a ratio of
 * about 30 for the rejection sampling. It must also stay small enough that the
 * extracted ||(d - d') * sigma_i|| <= 2 * (2 * sqrt(DEGREE) * SIGMA_S) = 32768
 * remains below sqrt(MODP) ~ 62501, the bound under which Lemma 1 guarantees
 * invertibility. */
#define SIGMA_S 256
