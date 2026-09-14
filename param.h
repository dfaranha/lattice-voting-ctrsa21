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
 * irreducible factors as Lemma 1 requires for delta = 2. Sized so that the
 * MSIS instance the binding property reduces to stays at 128 bits; see
 * PARAMS.md. Changing it requires new CRT constants P0, P1 in commit.c. */
#define MODP 	2553802997
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

/* Module dimensions. Table 1 of the paper ties the encryption dimensions to
 * the commitment ones, so they are all derived here rather than declared
 * separately in commit.h and encrypt.h, where they could drift apart. */

/* Width k of the commitment matrix. The hiding property is MLWE of rank
 * k - n - 1, so k = 3 leaves rank 1, which is only 72 bits; k = 4 gives
 * rank 2. There is nothing in between. */
#define WIDTH 	4
/* Height n of the commitment matrix. */
#define HEIGHT 	1
/* Maximum l1-norm nu of a challenge in C. The binding bound of the commitment
 * scheme is 16 * SIGMA_C * sqrt(NONZERO * DEGREE) and SIGMA_C itself grows
 * linearly in NONZERO, so the MSIS bound grows as NONZERO^1.5: keeping nu no
 * larger than the soundness error needs is worth more than it looks. At 20 the
 * challenge space is C(1024, 20) ~ 2^139, which covers the 4 * tau * t term of
 * Theorem 1 with room to spare. */
#define NONZERO 20
/* The infinity-norm bound beta_oo on the commitment randomness. */
#define BETA 	1
/* Dimension l of the encryption matrix, equal to k - n. */
#define DIM 	(WIDTH - HEIGHT)
/* Dimension kappa of the encryption message space, equal to the length k of
 * the commitment randomness that gets verifiably encrypted. */
#define VECTOR 	WIDTH

/* Standard deviation for discrete Gaussians masking the commitment randomness,
 * following the code's own convention SIGMA_C = 27.06 * NONZERO * sqrt(WIDTH *
 * DEGREE), which is about 1.23 times the bound in Table 1. */
#define SIGMA_C 34641
/* Standard deviation for the Gaussian masking the committed permutation
 * elements sigma_i of Lemma 5. These are monomials, so ||d * sigma_i|| equals
 * ||d|| <= sqrt(2 * NONZERO) ~ 6.3 exactly, and SIGMA_S = 190 leaves a ratio of
 * about 30 for the rejection sampling. It must also stay small enough that the
 * extracted ||(d - d') * sigma_i|| <= 2 * (2 * sqrt(DEGREE) * SIGMA_S) = 24320
 * remains below sqrt(MODP / 2) ~ 35735, the bound under which Lemma 1
 * guarantees invertibility. */
#define SIGMA_S 190
/* Standard deviation of the discrete Gaussian in the verifiable encryption.
 * Table 1 puts this at 11 * NONZERO * sqrt(VECTOR * DEGREE * (3 + BETA)) =
 * 28160, so the value below is conservative; it is kept unchanged because the
 * decryption bound of vericrypt has not been re-derived for a smaller one. */
#define SIGMA_E	54000

#endif /* PARAM_H */
