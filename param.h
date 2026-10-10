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
 * MSIS instance the binding property reduces to clears 128 bits even at the
 * conservative bound 2 ||delta||_1 * 4 sqrt(N) SIGMA_C of SOUNDNESS.md, Section
 * 5: 134.3 there, and 168.2 at the paper's 16 SIGMA_C sqrt(nu N); see
 * PARAMS.md. It was 2553802997, about 2^31.25, which gave 105.4 at the
 * conservative bound. Changing it requires new CRT constants P0, P1 in
 * commit.c, and a q in encrypt.c large enough for decryption. */
#define MODP 	137438953349
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
/* Parameter nu in the commitment scheme: the number of non-zero coefficients,
 * each +-1, of a challenge in C, and so its l1-norm. A difference of two
 * challenges, which is what an opening carries, has l1-norm at most
 * 2 * NONZERO. The binding bound of the commitment scheme is
 * 16 * SIGMA_C * sqrt(NONZERO * DEGREE) and SIGMA_C itself grows linearly in
 * NONZERO, so the MSIS bound grows as NONZERO^1.5: keeping nu no larger than
 * the soundness error needs is worth more than it looks. At 20 the challenge
 * space is 2^20 * C(1024, 20) ~ 2^158.7, which covers the 4 * tau * t term of
 * Theorem 1 with room to spare. */
#define NONZERO 20
/* The infinity-norm bound beta_oo on the commitment randomness. */
#define BETA 	1
/* Dimension l of the encryption matrix, equal to k - n. */
#define DIM 	(WIDTH - HEIGHT)
/* Dimension kappa of the encryption message space, equal to the length k of
 * the commitment randomness that gets verifiably encrypted. */
#define VECTOR 	WIDTH

/* Standard deviation for discrete Gaussians masking the commitment randomness.
 * A challenge d has NONZERO coefficients of +-1, so for ternary randomness r
 * the shift is ||d * r|| <= T = NONZERO * BETA * sqrt(WIDTH * DEGREE) = 1280.
 * Rejection sampling at the code's M = 1.75 keeps the statistical distance at
 * 2^-100 when alpha = SIGMA_C / T satisfies M >= exp(12 / alpha +
 * 1 / (2 * alpha^2)), that is alpha >= 21.485. So SIGMA_C = ceil(21.485 *
 * NONZERO * sqrt(WIDTH * DEGREE)), the smallest width that keeps the bound,
 * and alpha = 21.485. The constant takes the place of 11 in Table 1. */
#define SIGMA_C 27501
/* Standard deviation for the Gaussian masking the committed permutation
 * elements sigma_i of Lemma 5. These are monomials, so ||d * sigma_i|| equals
 * ||d|| = sqrt(NONZERO) ~ 4.47 exactly for a challenge d with NONZERO
 * coefficients of +-1, and SIGMA_S = ceil(21.485 * sqrt(NONZERO)) = 97 is the
 * smallest width that keeps the same 2^-100 bound at M = 1.75, at a ratio of
 * about 21.7. It must also stay small enough that the
 * extracted ||(d - d') * sigma_i|| <= 2 * (2 * sqrt(DEGREE) * SIGMA_S) = 12416
 * remains below sqrt(MODP / 2) ~ 262144, the bound under which Lemma 1
 * guarantees invertibility. */
#define SIGMA_S 97
/* Standard deviation of the discrete Gaussian in the verifiable encryption.
 * Table 1 puts this at 11 * NONZERO * sqrt(VECTOR * DEGREE * (3 + BETA)) =
 * 28160, so the value below is conservative; it is kept unchanged because the
 * decryption bound of vericrypt has not been re-derived for a smaller one. */
#define SIGMA_E	54000

#endif /* PARAM_H */
