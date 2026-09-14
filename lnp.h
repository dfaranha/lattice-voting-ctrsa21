/**
 * @defgroup lnp Lyubashevsky-Nguyen-Plancon proof machinery.
 */
/**
 * @file
 *
 * Interface of the LNP-style proof machinery needed by the is_bin sub-proof
 * of Protocol 1 in ePrint 2025/658. See ISBIN-PLAN.md for what this is for
 * and LNP-PARAMS.md for which parts are and are not established.
 *
 * @ingroup lnp
 */

#ifndef LNP_H
#define LNP_H

#include <flint/flint.h>
#include <flint/nmod_poly.h>

#include "param.h"
#include "commit.h"

/*============================================================================*/
/* Constant definitions                                                       */
/*============================================================================*/

/* Number of independent masks used to prove that a constant coefficient is
 * zero. Each one contributes a factor 1/MODP to the soundness error, so four
 * of them give about 2^-127 at this modulus. */
#define LNP_LAMBDA 	4

/* Number of message slots in the multi-slot commitment. Slots 0 to 3 are the
 * is_bin witness, its claimed product and two garbage terms, slot 4 is the
 * projection mask of the range proof, and the remaining LNP_LAMBDA slots hold
 * the masks of the constant-coefficient proof, so that one commitment serves
 * every part of the proof. */
#define SLOTS 	(5 + LNP_LAMBDA)

/* Number of coordinates the approximate range proof projects onto. The
 * projection lemma needs 256 of them for a 2^-128 soundness error. */
#define PROJ 	256

/* Ratio between the projection mask and the value it hides. Lower means a
 * tighter certified bound but more rejection-sampling repetitions; 9 gives
 * about 3.8 and leaves the certified bound a factor 1.78 inside the ceiling
 * that no-wraparound imposes. */
#define TAU_PROJ 	9

/* Slot holding the witness. */
#define SLOT_S 	0

/* Slot holding the value whose constant coefficient is proven to be zero. */
#define SLOT_F 	1

/* Slot holding the projection mask, packed into the first PROJ coefficients
 * of a ring element. */
#define SLOT_W 	4

/* First slot holding a constant-coefficient mask. */
#define SLOT_G 	5

/* Rank of the MLWE instance that hides the commitment, which is the number of
 * randomness components beyond those consumed by the Ajtai part and by the
 * message slots. Rank 1 is only about 72 bits; see LNP-PARAMS.md. */
#define LNP_RANK 	2

/* Width of the LNP commitment randomness. Unlike commit.h, which carries a
 * single message slot, a BDLOP key with SLOTS message rows in Hermite normal
 * form needs one randomness component per row on top of the Ajtai part. */
#define LNP_WIDTH 	(HEIGHT + SLOTS + LNP_RANK)

/*============================================================================*/
/* Type definitions                                                           */
/*============================================================================*/

/* A commitment key with SLOTS message rows instead of the single row that
 * commit.h provides. The Ajtai part B1 is shared. */
typedef struct _lnpkey_t {
	pcrt_poly_t B1[HEIGHT][LNP_WIDTH];
	pcrt_poly_t b2[SLOTS][LNP_WIDTH];
} lnpkey_t;

/* A commitment to SLOTS messages under one randomness vector. */
typedef struct _lnpcom_t {
	pcrt_poly_t c1[HEIGHT];
	pcrt_poly_t c2[SLOTS];
} lnpcom_t;

/* A proof that the committed messages satisfy a quadratic relation. The
 * garbage slots of the commitment are part of the proof rather than of the
 * statement, because they depend on the masking. */
typedef struct _lnpproof_t {
	pcrt_poly_t w[HEIGHT];			/* Ajtai part of the first message. */
	pcrt_poly_t t;					/* The masked challenge-free term. */
	pcrt_poly_t z[LNP_WIDTH];		/* The masked opening. */
} lnpproof_t;


/* The whole is_bin argument in one proof: the product relation, the constant
 * coefficient and the norm bound. All three speak about one commitment under
 * one randomness, so they share a single mask, a single challenge and a single
 * masked opening z, instead of carrying one each. */
typedef struct _lnpbinproof_t {
	pcrt_poly_t w[HEIGHT];			/* Ajtai part of the first message. */
	pcrt_poly_t t;					/* Product relation's masked term. */
	ulong zp[PROJ];					/* The masked projection. */
	pcrt_poly_t h[LNP_LAMBDA];		/* The aggregated values. */
	pcrt_poly_t v[LNP_LAMBDA];		/* The masked openings of the relations. */
	pcrt_poly_t z[LNP_WIDTH];		/* The single masked opening. */
} lnpbinproof_t;

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Apply the automorphism sigma_k : X -> X^k to a polynomial given in
 * coefficient representation. This is the reference implementation that the
 * CRT-aware one below is tested against.
 *
 * @param[out] c			- the result, may alias a.
 * @param[in] a				- the polynomial to map.
 * @param[in] k				- the exponent, which must be odd.
 */
void lnp_auto(nmod_poly_t c, const nmod_poly_t a, slong k);

/**
 * Apply the automorphism sigma_k to a polynomial in CRT representation.
 *
 * sigma_k maps the ideal generated by one irreducible factor onto the ideal
 * generated by the other when k = 3 mod 4, and onto itself when k = 1 mod 4,
 * so this is not a componentwise operation.
 *
 * @param[out] c			- the result, may alias a.
 * @param[in] a				- the polynomial in CRT representation.
 * @param[in] k				- the exponent, which must be odd.
 */
void lnp_auto_crt(pcrt_poly_t c, pcrt_poly_t a, slong k);

/**
 * Return whether sigma_k exchanges the two CRT components.
 *
 * @param[in] k				- the exponent, which must be odd.
 * @return 1 if the components are exchanged, 0 if they are fixed.
 */
int lnp_auto_swaps(slong k);

/**
 * Initialise, generate and free a multi-slot commitment key.
 */
void lnp_keyinit(lnpkey_t *key);
void lnp_keygen(lnpkey_t *key, flint_rand_t rand);
void lnp_keyfree(lnpkey_t *key);

/**
 * Initialise and free a commitment and a proof.
 */
void lnp_com_init(lnpcom_t *com);
void lnp_com_free(lnpcom_t *com);
void lnp_proof_init(lnpproof_t *pi);
void lnp_proof_free(lnpproof_t *pi);

/**
 * Commit to SLOTS messages under one randomness vector.
 *
 * @param[out] com			- the resulting commitment.
 * @param[in] m				- the messages, in CRT representation.
 * @param[in] key			- the commitment key.
 * @param[in] r				- the randomness, in CRT representation.
 */
void lnp_commit(lnpcom_t *com, pcrt_poly_t m[SLOTS], lnpkey_t *key,
		pcrt_poly_t r[LNP_WIDTH]);

/**
 * Prove that the committed messages satisfy m[0] * m[1] = m[2].
 *
 * Slot 3 carries the garbage term, so the caller supplies its randomness but
 * not its message: the prover computes it. The commitment passed in must
 * already hold slots 0 to 2.
 *
 * @param[out] pi			- the resulting proof.
 * @param[in,out] com		- the commitment, whose garbage slot is filled in.
 * @param[in] m				- the three messages, in CRT representation.
 * @param[in] key			- the commitment key.
 * @param[in] r				- the commitment randomness, in CRT representation.
 */
void lnp_quad_prover(lnpproof_t *pi, lnpcom_t *com, pcrt_poly_t m[3],
		lnpkey_t *key, pcrt_poly_t r[LNP_WIDTH]);

/**
 * Verify a proof that the committed messages satisfy m[0] * m[1] = m[2].
 *
 * @param[in] pi			- the proof.
 * @param[in] com			- the commitment.
 * @param[in] key			- the commitment key.
 * @return 1 if the proof is accepted, 0 otherwise.
 */
int lnp_quad_verifier(lnpproof_t *pi, lnpcom_t *com, lnpkey_t *key);

/**
 * Return the polynomial 1 + X + ... + X^(DEGREE-1) in CRT representation.
 *
 * The binary constraint on the coefficients of s is <s, s - ones> = 0, so this
 * is the vector of all ones seen as a ring element.
 */
void lnp_ones(pcrt_poly_t out);

/**
 * Compute the value an honest prover commits in slot SLOT_F, namely
 * sigma_{-1}(s) * (s - ones), whose constant coefficient is the sum over j of
 * s_j (s_j - 1).
 *
 * @param[out] f			- the product, in CRT representation.
 * @param[in] s				- the witness, in CRT representation.
 */
void lnp_isbin_product(pcrt_poly_t f, pcrt_poly_t s);


/**
 * Sample a polynomial uniformly at random subject to its constant coefficient
 * being zero. These are the masks of the constant-coefficient proof.
 *
 * @param[out] g			- the sampled polynomial, in CRT representation.
 * @param[in] rand			- the source of randomness.
 */
void lnp_sample_ct_zero(pcrt_poly_t g, flint_rand_t rand);



/**
 * Sample the projection mask: a ring element whose first PROJ coefficients are
 * Gaussian of width TAU_PROJ * sqrt(PROJ * DEGREE / 2) and whose remaining
 * coefficients are zero.
 *
 * @param[out] w			- the mask, in CRT representation.
 * @param[out] raw			- the same mask in coefficient representation.
 */
void lnp_sample_proj_mask(pcrt_poly_t w, nmod_poly_t raw);


/**
 * Expose the aggregation scalars, for tests that play the part of a prover
 * trying to adapt its masks to them.
 */
void lnp_scalars_for_test(ulong nu[LNP_LAMBDA], lnpkey_t *key, lnpcom_t *com,
		ulong z[PROJ]);

/**
 * Initialise and free the combined is_bin proof.
 */
void lnp_binproof_init(lnpbinproof_t *pi);
void lnp_binproof_free(lnpbinproof_t *pi);

/**
 * Prove the whole is_bin argument in one proof: that slot SLOT_F holds
 * sigma_{-1}(s) * (s - ones) for the witness s in slot SLOT_S, that its
 * constant coefficient is zero, and that s is short.
 *
 * All three speak about one commitment under one randomness, so they share a
 * single mask, a single challenge and a single masked opening rather than
 * carrying one each.
 *
 * @param[out] pi			- the resulting proof.
 * @param[in,out] com		- the commitment, whose garbage slots are filled in.
 * @param[in] s				- the witness, in CRT representation.
 * @param[in] f				- the claimed product, in CRT representation.
 * @param[in] w				- the projection mask, in coefficient form.
 * @param[in] g				- the constant-coefficient masks.
 * @param[in] key			- the commitment key.
 * @param[in] r				- the commitment randomness.
 * @return 1 if a transcript was produced, 0 if the projection was rejected.
 */
int lnp_bin_prover(lnpbinproof_t *pi, lnpcom_t *com, pcrt_poly_t s,
		pcrt_poly_t f, nmod_poly_t w, pcrt_poly_t g[LNP_LAMBDA],
		lnpkey_t *key, pcrt_poly_t r[LNP_WIDTH]);

/**
 * Verify the proof produced by lnp_bin_prover.
 */
int lnp_bin_verifier(lnpbinproof_t *pi, lnpcom_t *com, lnpkey_t *key);


#endif /* LNP_H */
