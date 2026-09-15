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
 * of them give about 2^-160 at this modulus. The 2^-127 this comment used to
 * quote was for the modulus before B5 raised it. */
#define LNP_LAMBDA 	4

/* Number of message slots in the multi-slot commitment: the is_bin witness,
 * its claimed product, and the projection mask of the range proof. Neither the
 * constant-coefficient masks nor the quadratic proof's garbage terms are here.
 * One set of each covers the whole batch, so they live in commitments of their
 * own rather than being paid for once per message. */
#define SLOTS 	3

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
#define SLOT_W 	2

/* Rank of the MLWE instance that hides the commitment, which is the number of
 * randomness components beyond those consumed by the Ajtai part and by the
 * message slots. Rank 1 is only about 72 bits; see LNP-PARAMS.md. */
#define LNP_RANK 	2

/* Width of the LNP commitment randomness. Unlike commit.h, which carries a
 * single message slot, a BDLOP key with SLOTS message rows in Hermite normal
 * form needs one randomness component per row on top of the Ajtai part. */
#define LNP_WIDTH 	(HEIGHT + SLOTS + LNP_RANK)

/* Width of the randomness of the batch-mask commitment. */
#define MASK_WIDTH 	(HEIGHT + LNP_LAMBDA + LNP_RANK)

/* The quadratic proof needs two garbage terms to cancel the challenge-linear
 * part of its expansion. Batched, one pair serves every message, so they get a
 * commitment of their own with two message rows. It cannot share the mask
 * commitment: that one has to be fixed before the projection scalars are
 * derived, whereas these depend on the batching challenge, which cannot be
 * drawn until every message's commitment exists. */
#define GARB_WIDTH 	(HEIGHT + 2 + LNP_RANK)

/*============================================================================*/
/* Type definitions                                                           */
/*============================================================================*/

/* A commitment key with SLOTS message rows instead of the single row that
 * commit.h provides. The Ajtai part B1 is shared. */
typedef struct _lnpkey_t {
	pcrt_poly_t B1[HEIGHT][LNP_WIDTH];
	pcrt_poly_t b2[SLOTS][LNP_WIDTH];
} lnpkey_t;

/* The key and commitment holding the constant-coefficient masks. There is one
 * of each for the whole batch, which is the point: the scalar aggregation pays
 * for LNP_LAMBDA masks however many statements it covers, so paying for them
 * once per message was waste. */
typedef struct _lnpmaskkey_t {
	pcrt_poly_t B1[HEIGHT][MASK_WIDTH];
	pcrt_poly_t b2[LNP_LAMBDA][MASK_WIDTH];
} lnpmaskkey_t;

typedef struct _lnpmaskcom_t {
	pcrt_poly_t c1[HEIGHT];
	pcrt_poly_t c2[LNP_LAMBDA];
} lnpmaskcom_t;

/* The key and commitment holding the batched garbage terms. */
typedef struct _lnpgarbkey_t {
	pcrt_poly_t B1[HEIGHT][GARB_WIDTH];
	pcrt_poly_t b2[2][GARB_WIDTH];
} lnpgarbkey_t;

typedef struct _lnpgarbcom_t {
	pcrt_poly_t c1[HEIGHT];
	pcrt_poly_t c2[2];
} lnpgarbcom_t;

/* The aggregated values, one set for the whole batch rather than one per
 * message, and the Ajtai first message that opens the mask commitment. The
 * latter is what binds the prover to the masks it committed: without it the
 * rows below are LNP_LAMBDA equations in MASK_WIDTH unknowns, and a prover
 * could pick any h with zero constant coefficient and solve for z_mask. */
typedef struct _lnpbatch_t {
	pcrt_poly_t w[HEIGHT];			/* Ajtai first message of the masks. */
	pcrt_poly_t h[LNP_LAMBDA];
	pcrt_poly_t v[LNP_LAMBDA];
	pcrt_poly_t gw[HEIGHT];			/* Ajtai first message of the garbage. */
	pcrt_poly_t T;					/* The batched challenge-free term. */
	pcrt_poly_t acc;				/* Verifier-side accumulator for T. */
} lnpbatch_t;

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
	ulong zp[PROJ];					/* The masked projection. */
} lnpbinproof_t;

/* What the setup phase derives from the commitment and the projection. It is
 * public, so the verifier rebuilds it rather than receiving it. */
typedef struct _lnpbinctx_t {
	pcrt_poly_t P[LNP_LAMBDA];
	pcrt_poly_t M[LNP_LAMBDA];
	ulong Z[LNP_LAMBDA];
	ulong nu[LNP_LAMBDA];
} lnpbinctx_t;

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

/*
 * The standalone quadratic proof that used to live here is gone. It needed
 * four slots, three for its messages and one for its garbage term, and SLOTS
 * is 3 now that the garbage terms are batched. The technique it demonstrated
 * is exercised in its real form by the is_bin tests, which drive the same
 * expansion through lnp_bin_first and lnp_bin_check.
 */


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
		lnpmaskcom_t *mcom, ulong z[PROJ]);

/**
 * Initialise and free the combined is_bin proof.
 */
void lnp_binproof_init(lnpbinproof_t *pi);
void lnp_binproof_free(lnpbinproof_t *pi);

/*
 * The is_bin argument, batched. One challenge and one opening are shared with
 * the linear proof of the shuffle, and one set of LNP_LAMBDA masks covers
 * every message, so the aggregated values h and v are accumulated across the
 * batch rather than produced per message.
 *
 * The order is: commit the masks, then per message setup and first, then the
 * caller's challenge, then per message check, then the batch check.
 */

/**
 * Initialise, generate and free the batch-mask key and commitment.
 */
void lnp_maskkey_init(lnpmaskkey_t *key);
void lnp_maskkey_gen(lnpmaskkey_t *key, flint_rand_t rand);
void lnp_maskkey_free(lnpmaskkey_t *key);
void lnp_maskcom_init(lnpmaskcom_t *com);
void lnp_maskcom_free(lnpmaskcom_t *com);
void lnp_mask_commit(lnpmaskcom_t *com, pcrt_poly_t g[LNP_LAMBDA],
		lnpmaskkey_t *key, pcrt_poly_t r[MASK_WIDTH]);

void lnp_batch_init(lnpbatch_t *b);
void lnp_batch_free(lnpbatch_t *b);
void lnp_batch_zero(lnpbatch_t *b);

/**
 * Initialise, generate and free the garbage-term key and commitment.
 */
void lnp_garbkey_init(lnpgarbkey_t *key);
void lnp_garbkey_gen(lnpgarbkey_t *key, flint_rand_t rand);
void lnp_garbkey_free(lnpgarbkey_t *key);
void lnp_garbcom_init(lnpgarbcom_t *com);
void lnp_garbcom_free(lnpgarbcom_t *com);
void lnp_garb_commit(lnpgarbcom_t *com, pcrt_poly_t g[2], lnpgarbkey_t *key,
		pcrt_poly_t r[GARB_WIDTH]);

/**
 * The garbage commitment's Ajtai first message, and the mask terms of the
 * batched challenge-free term T.
 */
void lnp_garb_first(lnpbatch_t *batch, lnpgarbkey_t *gkey,
		pcrt_poly_t yg[GARB_WIDTH]);

/**
 * Everything for one message that does not depend on the mask, accumulating
 * this message's share of the aggregated values h.
 *
 * @return 1 if a transcript can be produced, 0 if the projection was rejected.
 */
int lnp_bin_setup(lnpbinproof_t *pi, lnpbinctx_t *ctx, lnpbatch_t *batch,
		lnpcom_t *com, lnpmaskcom_t *mcom, pcrt_poly_t s, pcrt_poly_t f,
		nmod_poly_t w, lnpkey_t *key);

/**
 * Rebuild the public part of the setup, for the verifier.
 */
void lnp_bin_public(lnpbinctx_t *ctx, lnpbinproof_t *pi, lnpcom_t *com,
		lnpmaskcom_t *mcom, lnpkey_t *key);

/**
 * This message's mask-dependent first messages, accumulating its share of v,
 * its rho-weighted share of the two batched garbage terms, and its share of
 * the batched challenge-free term T.
 *
 * @param[in] rho			- this message's batching challenge.
 * @param[in,out] garb		- the two garbage accumulators for the batch.
 */
void lnp_bin_first(lnpbinctx_t *ctx, lnpbatch_t *batch, lnpcom_t *com,
		pcrt_poly_t s, lnpkey_t *key, pcrt_poly_t r[LNP_WIDTH],
		pcrt_poly_t y[LNP_WIDTH], pcrt_poly_t rho, pcrt_poly_t garb[2]);

/**
 * Accumulate this message's share of the two aggregated relations: the
 * rho-weighted quadratic terms into batch->acc, and the range and
 * constant-coefficient rows into acc.
 */
int lnp_bin_check(lnpbinproof_t *pi, lnpbinctx_t *ctx, lnpbatch_t *batch,
		lnpcom_t *com, lnpkey_t *key, pcrt_poly_t d, pcrt_poly_t z[LNP_WIDTH],
		pcrt_poly_t acc[LNP_LAMBDA], pcrt_poly_t rho);

/**
 * The batch-wide check, once every message has contributed.
 */
int lnp_batch_check(lnpbatch_t *batch, lnpmaskcom_t *mcom, lnpmaskkey_t *mkey,
		lnpgarbcom_t *gcom, lnpgarbkey_t *gkey, pcrt_poly_t d,
		pcrt_poly_t zm[MASK_WIDTH], pcrt_poly_t zg[GARB_WIDTH],
		pcrt_poly_t acc[LNP_LAMBDA]);

/**
 * The batch-mask commitment's contribution to the first message.
 */
void lnp_batch_first(lnpbatch_t *batch, lnpmaskkey_t *mkey,
		pcrt_poly_t ym[MASK_WIDTH]);

void lnp_binctx_init(lnpbinctx_t *ctx);
void lnp_binctx_free(lnpbinctx_t *ctx);


#endif /* LNP_H */
