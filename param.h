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

/* Modulus p defining the cyclotomic ring for the commitment scheme. Prime and
 * congruent to 5 mod 8, so that (x^DEGREE + 1) splits into exactly NCRT = 2
 * factors as Lemma 1 requires. Sized so that the approximate range proof can
 * certify a bound under sqrt(p/2); see LNP-PARAMS.md. Changing it requires new
 * CRT constants P0 and P1 in commit.c, and a new Q in encrypt.c. Changing
 * DEGREE does not: P0 is a square root of -1, which does not depend on the
 * degree, and x^(DEGREE/2) +- P0 stays irreducible. */
#define MODP 	1099511627917
/* Degree of the polynomial defining the cyclotomic ring for the commitment scheme. */
#ifndef DEGREE
#define DEGREE 	2048
#endif
/* Degree of each polynomial used to define the CRT representation. */
#define DEGCRT 	(DEGREE >> 1)
/* Number of irreducible factors of (x^DEGREE + 1), and therefore the number of
 * components of a CRT representation. This is a property of the ring, and is
 * unrelated to the module dimensions DIM, WIDTH and HEIGHT: it only happens to
 * share the value 2 with DIM, which made several loops index a CRT component
 * with a module bound. */
#define NCRT 	2
/* Standard deviation for discrete Gaussians: about 23.44 * NONZERO *
 * sqrt(WIDTH * DEGREE), so it grows as sqrt(DEGREE). */
#define SIGMA_C 76368
/* Standard deviation for the Gaussian masking the committed permutation
 * elements sigma_i of Lemma 5. These are monomials, so ||d * sigma_i|| equals
 * ||d|| <= sqrt(2 * NONZERO) ~ 8.5 exactly regardless of DEGREE, which is why
 * this one does not scale with it. SIGMA_S = 256 leaves a ratio of about 30
 * for the rejection sampling. It must also stay small enough that the
 * extracted ||(d - d') * sigma_i|| <= 2 * (2 * sqrt(DEGREE) * SIGMA_S), which
 * is 46341 at DEGREE 2048, remains below sqrt(MODP) ~ 1048576, the bound under
 * which Lemma 1 guarantees invertibility. */
#define SIGMA_S 256
/* Standard deviation for the masks of the batched proof. One challenge shared
 * across MSGS statements means one rejection test over their concatenation, so
 * the term being masked is sqrt(MSGS) times longer than for a single statement
 * and the mask has to widen to match. */
#define SIGMA_B 381840
/* Standard deviation of the projection mask in the approximate range proof:
 * TAU_PROJ * sqrt(PROJ * DEGREE / 2), the width needed to hide a projection of
 * an honest binary witness. The bound it certifies grows with it, so the
 * margin against the wraparound ceiling sqrt(p/2) narrows as DEGREE rises:
 * 1.78x at DEGREE 1024 and 1.26x at 2048. */
#define SIGMA_P 4608
// Standard deviation of discrete Gaussian
#define SIGMA_E	54000

#endif /* PARAM_H */
