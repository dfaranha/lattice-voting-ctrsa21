/**
 * @defgroup encrypt Lattice-based encryption.
 */
/**
 * @file
 *
 * Interface of the lattice-based encryption scheme.
 *
 * @ingroup encrypt
 */

#include <flint/flint.h>
#include <flint/fmpz_mod_poly.h>

/*============================================================================*/
/* Constant definitions                                                       */
/*============================================================================*/

/* The dimension of the encryption matrix. */
#define DIM		2
/* The dimension of the linear relation matrix. */
#define LAMBDA	1

/*============================================================================*/
/* Type definitions                                                           */
/*============================================================================*/

/* Type that represents a polynomial in CRT representation. */
typedef fmpz_mod_poly_t qcrt_poly_t[2];

/* Type that represents a public key for the encryption scheme. */
typedef struct _publickey_t {
	qcrt_poly_t A[DIM][DIM];
	qcrt_poly_t t[DIM];
} publickey_t;

/* Type that represents a private key for the encryption scheme. */
typedef struct _privatekey_t {
	qcrt_poly_t s1[DIM];
	qcrt_poly_t s2[DIM];
} privatekey_t;

/* Type that represents a ciphertext for the encryption scheme. */
typedef struct _ciphertext_t {
	qcrt_poly_t v[DIM];
	qcrt_poly_t w;
} ciphertext_t;

/* Recover polynomial from CRT representation.
 *
 * @param[in] c 		- the resulting polynomial.
 * @param[in] a 		- the polynomial in CRT representation.
 */
void qcrt_poly_rec(fmpz_mod_poly_t c, qcrt_poly_t a);

/**
 * Initialize the commitment module.
 */
void encrypt_setup();

/**
 * Finalize the commitment module.
 */
void encrypt_finish();

/**
 * Sample a short polynomial.
 *
 * @param[out] r		- the polynomial to sample.
 * @param[in] ctx 		- the context for modular arithmetic.
 */
void encrypt_sample_short(fmpz_mod_poly_t r, fmpz_mod_ctx_t ctx);

/**
 * Sample a short polynomial in CRT representation.
 *
 * @param[out] r		- the polynomial to sample.
 * @param[in] ctx 		- the context for modular arithmetic.
 */
void encrypt_sample_short_crt(qcrt_poly_t r, fmpz_mod_ctx_t ctx);

/**
 * Returns the small modulus p for the encryption system.
 *
 * @return the small modulus.
 */
fmpz_t *encrypt_modulus();

/**
 * Returns the large modulus q for the encryption system.
 *
 * @return the large modulus.
 */
fmpz_t *encrypt_large_modulus();

/**
 * Returns the context for the small modulus p for the encryption system.
 *
 * @return the context for the small modulus.
 */
fmpz_mod_ctx_t *encrypt_modulus_ctx();

/**
 * Returns the context for the large modulus p for the encryption system.
 *
 * @return the context for the large modulus.
 */
fmpz_mod_ctx_t *encrypt_large_modulus_ctx();

/**
 * Return the i-th polynomial defining the CRT representation.
 * The returned polynomial is such that it factors the polynomial defining Rp.
 *
 * @return the i-th irreducible polynomial defining the CRT representation.
 */
fmpz_mod_poly_t *encrypt_irred(int i);

/**
 * Return the i-th polynomial defining the CRT representation.
 * The returned polynomial is such that it factors the polynomial defining Rq.
 *
 * @return the i-th irreducible polynomial defining the CRT representation.
 */
fmpz_mod_poly_t *encrypt_large_irred(int i);

/**
 * Return the polynomial (x^N + 1) defining the cyclotomic ring Rq.
 *
 * @return the polynomial defining the cyclotomic ring Rq.
 */
fmpz_mod_poly_t *encrypt_poly();

/**
 * Generate a key pair for the encryption scheme using a PRNG.
 *
 * @param[out] pk 		- the generated public key.
 * @param[out] sk 		- the generated private key.
 * @param[in] rand		- the random number generator.
 */
void encrypt_keygen(publickey_t *pk, privatekey_t *sk, flint_rand_t rand);

/**
 * Free a key pair for the encryption scheme.
 *
 * @param[in] pk 		- the public key to free.
 * @param[in] sk		- the private key to free.
 */
void encrypt_keyfree(publickey_t *pk, privatekey_t *sk);

/**
 * Encrypt a message under a public key.
 *
 * @param[out]			- the resulting ciphertext.
 * @param[in]			- the message to encrypt.
 * @param[in]			- the public key to encrypt under.
 * @param[in] rand		- the random number generator.
 */
void encrypt_doit(ciphertext_t *c, fmpz_mod_poly_t m, publickey_t *pk, flint_rand_t rand);

/**
 * Internal function to encrypt a message under a public key.
 *
 * @param[out]			- the resulting ciphertext.
 * @param[in]			- the message to encrypt.
 * @param[in]			- the public key to encrypt under.
 * @param[in] rand		- the random number generator.
 */
void encrypt_make(ciphertext_t *c, qcrt_poly_t r[DIM], qcrt_poly_t e[DIM],
	qcrt_poly_t e_, fmpz_mod_poly_t m, publickey_t *pk);

/**
 * Decrypt a ciphertext into the original plaintext message.
 *
 * @param[out] m		- the resulting message.
 * @param[in] chall		- an optional decryption challenge.
 * @param[in] c			- the ciphertext to decrypt.
 * @param[in] sk		- the private key to decrypt.
 */
int encrypt_undo(fmpz_mod_poly_t m, fmpz_mod_poly_t chall, ciphertext_t *c, privatekey_t *sk);

/**
 * Free a ciphertext.
 *
 * @param[in] c			- the ciphertext to free.
 */
void encrypt_free(ciphertext_t *c);
