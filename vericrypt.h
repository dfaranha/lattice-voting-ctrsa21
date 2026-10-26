/**
 * @file
 *
 * Interface of the verifiable encryption: a Lyubashevsky-Neven proof that a
 * ciphertext encrypts a short witness m of a linear relation <t, m> = u over
 * R_p, used by a voter to encrypt the randomness of their commitment.
 *
 * @ingroup encrypt
 */

#ifndef VERICRYPT_H
#define VERICRYPT_H

#include "param.h"
#include "encrypt.h"
#include "sha.h"

/* Type that represents a verifiable encryption ciphertext. */
typedef struct _veritext_t {
	ciphertext_t cipher[VECTOR];
	fmpz_mod_poly_t c;
	fmpz_mod_poly_t r[VECTOR][DIM][NCRT];
	fmpz_mod_poly_t e[VECTOR][DIM][NCRT];
	fmpz_mod_poly_t e_[VECTOR][NCRT];
	fmpz_mod_poly_t u[VECTOR];
} veritext_t;

/**
 * Initialise a verifiable ciphertext. Must be called before vericrypt_doit,
 * and released with vericrypt_free.
 *
 * @param[out] out		- the ciphertext to initialise.
 */
void vericrypt_init(veritext_t *out);

/**
 * Free a verifiable ciphertext.
 *
 * @param[out] out		- the ciphertext to free.
 */
void vericrypt_free(veritext_t *out);

/**
 * Encrypt the short witness m of <t, m> = u under pk, with a proof.
 *
 * @param[out] out		- the verifiable ciphertext.
 * @param[in] t			- the relation.
 * @param[in] u			- its right-hand side.
 * @param[in] m			- the witness, short, modulo p.
 * @param[in] pk		- the public key.
 * @param[in] rand		- FLINT's random state.
 * @return 1 on success, 0 otherwise.
 */
int vericrypt_doit(veritext_t *out, fmpz_mod_poly_t t[VECTOR],
		fmpz_mod_poly_t u, fmpz_mod_poly_t m[VECTOR], publickey_t *pk,
		flint_rand_t rand);

/**
 * Verify a verifiable ciphertext for the relation (t, u) under pk.
 *
 * @return 1 if it verifies, 0 otherwise.
 */
int vericrypt_verify(veritext_t *in, fmpz_mod_poly_t t[VECTOR],
		fmpz_mod_poly_t u, publickey_t *pk);

/**
 * Decrypt a verifiable ciphertext, after verifying it. What comes out is a
 * relaxed witness: m with <t, m> = c u for a difference of challenges c, so
 * that for an honest prover m is c times the encrypted witness.
 *
 * @param[out] m		- the relaxed witness.
 * @param[out] c		- the difference of challenges.
 * @return 1 on success, 0 if the ciphertext does not verify or decrypt.
 */
int vericrypt_undo(fmpz_mod_poly_t m[VECTOR], fmpz_mod_poly_t c,
		veritext_t *in, fmpz_mod_poly_t t[VECTOR], fmpz_mod_poly_t u,
		publickey_t *pk, privatekey_t *sk);

#endif /* VERICRYPT_H */
