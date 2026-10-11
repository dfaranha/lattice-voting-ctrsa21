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
 * certify a bound under sqrt(p/2); see LNP-PARAMS.md. It is the largest such
 * prime below 2^40, 2^40 - 195, so that p - 1 fits in 40 bits: a uniform
 * coefficient packs in 40 and hashes in 5 bytes. Changing it requires new
 * CRT constants P0 and P1 in commit.c, and a new Q in encrypt.c. Changing
 * DEGREE does not: P0 is a square root of -1, which does not depend on the
 * degree, and x^(DEGREE/2) +- P0 stays irreducible. */
#define MODP 	1099511627581
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
/* Standard deviation for discrete Gaussians masking one commitment opening
 * d * r, r ternary in WIDTH ring elements, so ||d * r|| <= T = NONZERO *
 * sqrt(WIDTH * DEGREE). SIGMA_C = ceil(21.485 * T): alpha = sigma / T = 21.485
 * is the least ratio at which rejection sampling with the code's M = 1.921
 * keeps its 2^-128 statistical bound, M = exp(14 / alpha + 1 / (2 alpha^2))
 * (LNP22, Lemma 2.14).
 * It grows as sqrt(DEGREE) and shrinks with WIDTH. On this branch nothing is
 * masked at this width any more; it sets the opening bound commit_open
 * accepts, and the unit SIGMA_B is derived from. */
#define SIGMA_C 30314
/* Standard deviation for the Gaussian masking the committed permutation
 * elements sigma_i of Lemma 5. These are monomials, so ||d * sigma_i|| equals
 * ||d|| = sqrt(NONZERO) ~ 4.24 exactly regardless of DEGREE, d having NONZERO
 * coefficients of +-1, which is why this one does not scale with it.
 * SIGMA_S = ceil(21.485 * sqrt(NONZERO)) = 92, a ratio of about 21.7, the
 * rigorous minimum for M = 1.921. It must also stay small enough that the
 * extracted ||(d - d') * sigma_i|| <= 2 * (2 * sqrt(DEGREE) * SIGMA_S), which
 * is 16654 at DEGREE 2048, remains below sqrt(MODP) ~ 1048576, the bound under
 * which Lemma 1 guarantees invertibility. The relaxed proof that used it was
 * replaced by is_bin, so no proof on this branch samples at this width. */
#define SIGMA_S 92
/* Standard deviation for the masks of the batched proof. One challenge shared
 * across MSGS statements means one rejection test over their concatenation:
 * per message the openings d * r and d * _r of the two WIDTH-wide commitments
 * and d * s of the LNP_WIDTH-wide one, and once for the batch those of the
 * mask and garbage commitments, MASK_WIDTH + GARB_WIDTH more. That is
 * K = MSGS * (2 WIDTH + LNP_WIDTH) + MASK_WIDTH + GARB_WIDTH = 285 ternary ring
 * elements, so T = NONZERO * sqrt(K * DEGREE) and the width has to be
 * 21.485 * T = SIGMA_C * sqrt(K / WIDTH), rounded up. shuffle.c asserts it. */
#define SIGMA_B 295465
/* Standard deviation of the projection mask in the approximate range proof:
 * TAU_PROJ * sqrt(PROJ * DEGREE / 2), the width needed to hide a projection of
 * an honest binary witness. The bound it certifies grows with it, so the
 * margin against the wraparound ceiling sqrt(p/2) narrows as DEGREE rises:
 * 1.78x at DEGREE 1024 and 1.26x at 2048. */
#define SIGMA_P 4608
// Standard deviation of discrete Gaussian
#define SIGMA_E	54000

#endif /* PARAM_H */
