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

/*============================================================================*/
/* Constant definitions                                                       */
/*============================================================================*/

/* Modulus p defining the cyclotomic ring for the commitment scheme. */
#define MODP 	3906450253
/* Degree of the polynomial defining the cyclotomic ring for the commitment scheme. */
#define DEGREE 	1024
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
// Standard deviation of discrete Gaussian
#define SIGMA_E	54000

#endif /* PARAM_H */
