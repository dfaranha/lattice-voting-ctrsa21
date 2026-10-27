/**
 * @file
 *
 * The voter's proof that its commitment holds a valid vote, on LaZer's linear
 * proofs: com = (c1, c2) = (sum_j B1_j r_j, sum_j b2_j r_j + m) with short r and
 * m = X^j or 0, bound to a context such as the election and the voter. The
 * verifiable encryption of r shows that the authority can open the
 * commitment; this proof shows that what it opens to is one vote.
 *
 * C, so that the C++ build of the rest can call it without compiling LaZer's
 * headers: elements of R_p are passed as their VALIDITY_DEG coefficients,
 * centred mod VALIDITY_P, as full polynomials rather than in CRT form.
 */

#ifndef VALIDITY_H
#define VALIDITY_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* The ring the commitments live in: DEGREE and MODP of param.h, which
 * voting.c checks. */
#define VALIDITY_DEG 	1024
#define VALIDITY_P 		549755813869LL
/* The commitment width, WIDTH of param.h. */
#define VALIDITY_WIDTH 	4

/* Initialise LaZer; call once before proving or verifying. */
void validity_setup(void);

/* An upper bound on the size in bytes of a proof; its encoding varies. */
size_t validity_proof_len(void);

/**
 * Prove that c1 = sum_j B1[j] r[j] and c2 = sum_j b2[j] r[j] + m with r ternary
 * and short and m = X^j or 0, bound to ctx. Refuses (returns 0, no proof) if
 * the witness is not of that form or the equations do not hold.
 *
 * @param[out] proof		- room for validity_proof_len() bytes.
 * @param[out] len			- the bytes written, at most that.
 * @return 1 if a proof was made, 0 otherwise.
 */
int validity_prove(uint8_t *proof, size_t *len, const int64_t *c1,
		const int64_t *c2, const int64_t *B1, const int64_t *b2,
		const int64_t *r, const int64_t *m, const uint8_t *ctx, size_t ctxlen);

/**
 * Verify a validity proof for the commitment (c1, c2) under the key (B1, b2)
 * and the context ctx. B1, b2 and r hold VALIDITY_WIDTH elements in turn.
 *
 * @return 1 if it verifies, 0 otherwise.
 */
int validity_verify(const uint8_t *proof, size_t len, const int64_t *c1,
		const int64_t *c2, const int64_t *B1, const int64_t *b2,
		const uint8_t *ctx, size_t ctxlen);

#ifdef __cplusplus
}
#endif

#endif /* VALIDITY_H */
