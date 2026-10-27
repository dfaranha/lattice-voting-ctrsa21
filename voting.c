/**
 * @file
 *
 * The voting protocol end to end, as the paper's application of the proof of
 * shuffle puts the pieces together, with MSGS voters:
 *
 *   1. Setup: a commitment key, and a key pair for the verifiable encryption,
 *      whose secret key the decryption authority holds.
 *   2. Voting: each voter commits to its vote, X^j for one of NCAND
 *      candidates, as com = (c1, c2) = (<B1, r>, <b2, r> + m), and verifiably
 *      encrypts the commitment randomness r: a proof that the ciphertexts hold
 *      a short r with <B1, r> = c1, the linear relation of the commitment.
 *      It also proves the vote valid (validity.c, on LaZer): the commitment is
 *      to m = X^j or 0 under short r, bound to the election and the voter.
 *   3. The board verifies every ballot and its validity proof, against its
 *      own commitment.
 *   4. The authority decrypts every ballot. Decryption gives a relaxed
 *      witness c r for a difference of challenges c, so it multiplies by
 *      c^-1 in R_p to recover r, checks the opening of the commitment, and
 *      reads the vote off it. It then shuffles the votes and publishes them
 *      with a proof of shuffle, as bytes.
 *   5. Anyone verifies the proof from those bytes, the commitments and the
 *      published votes, and tallies.
 *
 * Then the same pipeline is attacked: a ballot copied onto another voter's
 * commitment, a published vote changed, a vote substituted with an honest
 * proof, and a corrupted proof must all be rejected.
 *
 * Each step uses only what its party has: voters their own vote and
 * randomness, the board the public values, the authority the ballots and its
 * key, the verifier the published bytes. The voters' votes are kept apart, as
 * a test oracle only.
 *
 * The prover and verifier of the proof of shuffle are internal to shuffle.c,
 * so this file includes it, with SHUFFLE_NO_MAIN leaving out its own tests.
 */

#define SHUFFLE_NO_MAIN
#include "shuffle.c"

#include "encrypt.h"
#include "vericrypt.h"
#include "validity.h"
#include <time.h>

static_assert(VALIDITY_DEG == DEGREE && VALIDITY_P == MODP && VALIDITY_WIDTH == WIDTH,
		"the validity proof is for this ring and commitment width");

/* Wall-clock seconds, for the phase timings the election reports. */
static double now_s(void) {
	struct timespec t;

	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + t.tv_nsec / 1e9;
}

/* Number of candidates: a vote for candidate j is the monomial X^j. */
#ifndef NCAND
#define NCAND 	8
#endif

/* A uniform integer in [0, n), from the operating system. */
static uint32_t uniform_below(uint32_t n) {
	uint32_t x, lim = UINT32_MAX - UINT32_MAX % n;
	do {
		randombytes((unsigned char *) &x, sizeof(x));
	} while (x >= lim);
	return x % n;
}

/* An element of R_p, from CRT to the encryption's representation modulo p. */
static void crt_to_mod(fmpz_mod_poly_t out, pcrt_poly_t in) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	pcrt_poly_rec(t, in);
	fmpz_mod_poly_zero(out, *encrypt_modulus_ctx());
	for (int i = 0; i < DEGREE; i++) {
		fmpz_mod_poly_set_coeff_ui(out, i, nmod_poly_get_coeff_ui(t, i),
				*encrypt_modulus_ctx());
	}
	nmod_poly_clear(t);
}

/* An element of R_p, from the encryption's representation to CRT. */
static void mod_to_crt(pcrt_poly_t out, fmpz_mod_poly_t in) {
	nmod_poly_t t;
	fmpz_t c;

	nmod_poly_init(t, MODP);
	fmpz_init(c);
	for (int i = 0; i < DEGREE; i++) {
		fmpz_mod_poly_get_coeff_fmpz(c, in, i, *encrypt_modulus_ctx());
		nmod_poly_set_coeff_ui(t, i, fmpz_get_ui(c));
	}
	for (int k = 0; k < NCRT; k++) {
		pcrt_poly_reduce(out[k], t, k);
	}
	fmpz_clear(c);
	nmod_poly_clear(t);
}

/* The coefficients of a in R_p, centred, as the validity proof takes them. */
static void nmod_to_i64(int64_t out[DEGREE], const nmod_poly_t a) {
	for (int i = 0; i < DEGREE; i++) {
		uint64_t c = i < nmod_poly_length(a) ? nmod_poly_get_coeff_ui(a, i) : 0;
		out[i] = c > MODP / 2 ? (int64_t) c - MODP : (int64_t) c;
	}
}

/* The same for an element in CRT form. */
static void crt_to_i64(int64_t out[DEGREE], pcrt_poly_t a) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	pcrt_poly_rec(t, a);
	nmod_to_i64(out, t);
	nmod_poly_clear(t);
}

/* A commitment and its key as the validity proof takes them: c1, c2, then
 * B1[0][j] and b2[j] for every j. */
typedef struct {
	int64_t c1[DEGREE], c2[DEGREE], B1[WIDTH * DEGREE], b2[WIDTH * DEGREE];
} validity_stmt_t;

static void validity_stmt(validity_stmt_t *st, commit_t *com, commitkey_t *key) {
	crt_to_i64(st->c1, com->c1);
	crt_to_i64(st->c2, com->c2);
	for (int j = 0; j < WIDTH; j++) {
		crt_to_i64(st->B1 + j * DEGREE, key->B1[0][j]);
		crt_to_i64(st->b2 + j * DEGREE, key->b2[j]);
	}
}

/* The context a validity proof is bound to: the election and the voter. */
#define VCTX_BYTES 	(16 + 16 + 4)
static void validity_ctx(uint8_t ctx[VCTX_BYTES], const uint8_t eid[16],
		int voter) {
	memcpy(ctx, "ctrsa21/validity", 16);
	memcpy(ctx + 16, eid, 16);
	for (int b = 0; b < 4; b++) {
		ctx[32 + b] = (uint8_t) ((uint32_t) voter >> (8 * b));
	}
}

/* A ballot's validity proof: proof[0 .. len - 1]. */
typedef struct {
	uint8_t *proof;
	size_t len;
} vproof_t;

/* The relation a ballot proves, for a commitment: t = B1 and u = c1. */
static void ballot_relation(fmpz_mod_poly_t t[VECTOR], fmpz_mod_poly_t u,
		commit_t *com, commitkey_t *key) {
	for (int j = 0; j < VECTOR; j++) {
		crt_to_mod(t[j], key->B1[0][j]);
	}
	crt_to_mod(u, com->c1);
}

/* A voter: commit to the vote, verifiably encrypt the randomness, and prove
 * the vote valid, bound to ctx. */
static int cast(commit_t *com, veritext_t *ballot, vproof_t *vp, nmod_poly_t vote,
		commitkey_t *key, publickey_t *pk, const uint8_t ctx[VCTX_BYTES],
		flint_rand_t rand) {
	pcrt_poly_t r[WIDTH];
	fmpz_mod_poly_t t[VECTOR], u, w[VECTOR];
	static validity_stmt_t st;
	static int64_t r64[WIDTH * DEGREE], m64[DEGREE];
	int ok;

	for (int j = 0; j < WIDTH; j++) {
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(r[j][k], MODP);
		}
		commit_sample_short_crt(r[j]);
	}
	commit_doit(com, vote, key, r);
	fmpz_mod_poly_init(u, *encrypt_modulus_ctx());
	for (int j = 0; j < VECTOR; j++) {
		fmpz_mod_poly_init(t[j], *encrypt_modulus_ctx());
		fmpz_mod_poly_init(w[j], *encrypt_modulus_ctx());
		crt_to_mod(w[j], r[j]);
	}
	ballot_relation(t, u, com, key);
	ok = vericrypt_doit(ballot, t, u, w, pk, rand);
	validity_stmt(&st, com, key);
	for (int j = 0; j < WIDTH; j++) {
		crt_to_i64(r64 + j * DEGREE, r[j]);
	}
	nmod_to_i64(m64, vote);
	ok &= validity_prove(vp->proof, &vp->len, st.c1, st.c2, st.B1, st.b2, r64,
			m64, ctx, VCTX_BYTES);
	for (int j = 0; j < VECTOR; j++) {
		fmpz_mod_poly_clear(t[j], *encrypt_modulus_ctx());
		fmpz_mod_poly_clear(w[j], *encrypt_modulus_ctx());
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(r[j][k]);
		}
	}
	fmpz_mod_poly_clear(u, *encrypt_modulus_ctx());
	return ok;
}

/* The validity proof of a commitment, under the context ctx. */
static int check_validity(vproof_t *vp, commit_t *com, commitkey_t *key,
		const uint8_t ctx[VCTX_BYTES]) {
	static validity_stmt_t st;

	validity_stmt(&st, com, key);
	return validity_verify(vp->proof, vp->len, st.c1, st.c2, st.B1, st.b2, ctx,
			VCTX_BYTES);
}

/* The board: a ballot and its validity proof verify against the commitment it
 * was cast with. */
static int check_ballot(veritext_t *ballot, vproof_t *vp, commit_t *com,
		commitkey_t *key, publickey_t *pk, const uint8_t ctx[VCTX_BYTES]) {
	fmpz_mod_poly_t t[VECTOR], u;
	int ok;

	fmpz_mod_poly_init(u, *encrypt_modulus_ctx());
	for (int j = 0; j < VECTOR; j++) {
		fmpz_mod_poly_init(t[j], *encrypt_modulus_ctx());
	}
	ballot_relation(t, u, com, key);
	ok = vericrypt_verify(ballot, t, u, pk) && check_validity(vp, com, key, ctx);
	for (int j = 0; j < VECTOR; j++) {
		fmpz_mod_poly_clear(t[j], *encrypt_modulus_ctx());
	}
	fmpz_mod_poly_clear(u, *encrypt_modulus_ctx());
	return ok;
}

/* The authority: decrypt a ballot, recover the exact randomness r = c^-1 (c r)
 * of the commitment, check the opening, and read the vote m = c2 - <b2, r>.
 * Returns 0 if any step fails. */
static int open_ballot(nmod_poly_t m, nmod_poly_t r[WIDTH][2],
		veritext_t *ballot, commit_t *com, commitkey_t *key, publickey_t *pk,
		privatekey_t *sk) {
	fmpz_mod_poly_t t[VECTOR], u, w[VECTOR], c;
	pcrt_poly_t cc, inv, one, wc, acc;
	nmod_poly_t tmp;
	int ok;

	fmpz_mod_poly_init(u, *encrypt_modulus_ctx());
	fmpz_mod_poly_init(c, *encrypt_modulus_ctx());
	for (int j = 0; j < VECTOR; j++) {
		fmpz_mod_poly_init(t[j], *encrypt_modulus_ctx());
		fmpz_mod_poly_init(w[j], *encrypt_modulus_ctx());
	}
	nmod_poly_init(tmp, MODP);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_init(cc[k], MODP);
		nmod_poly_init(inv[k], MODP);
		nmod_poly_init(one[k], MODP);
		nmod_poly_init(wc[k], MODP);
		nmod_poly_init(acc[k], MODP);
		nmod_poly_one(one[k]);
	}

	ballot_relation(t, u, com, key);
	ok = vericrypt_undo(w, c, ballot, t, u, pk, sk);
	if (ok) {
		/* c is a short difference of challenges, so it is invertible in R_p
		 * (Lemma 1), which means invertible in both CRT components. */
		mod_to_crt(cc, c);
		for (int k = 0; k < NCRT && ok; k++) {
			ok = nmod_poly_invmod(inv[k], cc[k], *commit_irred(k));
		}
	}
	for (int j = 0; j < WIDTH && ok; j++) {
		mod_to_crt(wc, w[j]);
		for (int k = 0; k < NCRT; k++) {
			pcrt_poly_mulmod(r[j][k], wc[k], inv[k], k);
		}
	}
	if (ok) {
		/* m = c2 - <b2, r>, and the opening (m, r) must check: r short and
		 * both rows of the commitment matching. */
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_zero(acc[k]);
			for (int j = 0; j < WIDTH; j++) {
				pcrt_poly_mulmod(tmp, key->b2[j][k], r[j][k], k);
				nmod_poly_add(acc[k], acc[k], tmp);
			}
			nmod_poly_sub(acc[k], com->c2[k], acc[k]);
		}
		pcrt_poly_rec(m, acc);
		ok = commit_open(com, m, key, r, one);
	}

	for (int j = 0; j < VECTOR; j++) {
		fmpz_mod_poly_clear(t[j], *encrypt_modulus_ctx());
		fmpz_mod_poly_clear(w[j], *encrypt_modulus_ctx());
	}
	fmpz_mod_poly_clear(u, *encrypt_modulus_ctx());
	fmpz_mod_poly_clear(c, *encrypt_modulus_ctx());
	nmod_poly_clear(tmp);
	for (int k = 0; k < NCRT; k++) {
		nmod_poly_clear(cc[k]);
		nmod_poly_clear(inv[k]);
		nmod_poly_clear(one[k]);
		nmod_poly_clear(wc[k]);
		nmod_poly_clear(acc[k]);
	}
	return ok;
}

/* The candidate a vote is for, or -1 if it is not X^j for some j < NCAND. */
static int candidate(nmod_poly_t m) {
	int j = -1;

	for (int i = 0; i < DEGREE; i++) {
		ulong c = nmod_poly_get_coeff_ui(m, i);
		if (c == 1 && i < NCAND && j < 0) {
			j = i;
		} else if (c != 0) {
			return -1;
		}
	}
	return j;
}

/* The proof as published: the bytes proof_walk writes, which is everything the
 * prover sends; the statement (the commitments and the published votes) and
 * the key are public. About 51 KB per message at these parameters, so it is
 * sized by MSGS: a fixed 8 MB held only about 160 messages, and prove() then
 * reported the overflow as an empty proof. */
#define PUBLISHED_CAP 	((size_t) MSGS * (64 << 10) + (1 << 20))
static uint8_t published[PUBLISHED_CAP];

/* The authority: prove that _m is a shuffle of the votes in com, given their
 * openings (m, r) and the permutation elements sigma_i = g(pi(i)), and write
 * the proof out. Returns its size in bytes, or 0 if it did not fit. */
static size_t prove(commit_t com[MSGS], nmod_poly_t m[MSGS],
		nmod_poly_t _m[MSGS], nmod_poly_t sigma[MSGS],
		nmod_poly_t r[MSGS][WIDTH][2], commitkey_t *key) {
	static commit_t d[MSGS], p[MSGS];
	static nmod_poly_t s[MSGS], u[MSGS][2];
	static nmod_poly_t y[MSGS][WIDTH][2], w[MSGS][WIDTH][2], _y[MSGS][WIDTH][2];
	static nmod_poly_t ys[MSGS][1][2], vs[MSGS][2];
	static nmod_poly_t t[MSGS][2], tp[MSGS][2], _t[MSGS][2];
	static uint8_t digest[MSGS][SHA256HashSize];
	bitwriter_t bw;
	size_t len;

	for (int i = 0; i < MSGS; i++) {
		commit_init(&d[i]);
		commit_init(&p[i]);
		nmod_poly_init(s[i], MODP);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(t[i][k], MODP);
			nmod_poly_init(tp[i][k], MODP);
			nmod_poly_init(_t[i][k], MODP);
			nmod_poly_init(ys[i][0][k], MODP);
			nmod_poly_init(vs[i][k], MODP);
			nmod_poly_init(u[i][k], MODP);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_init(y[i][j][k], MODP);
				nmod_poly_init(w[i][j][k], MODP);
				nmod_poly_init(_y[i][j][k], MODP);
			}
		}
	}

	shuffle_prover(y, w, _y, ys, t, tp, _t, vs, u, d, s, com, p, m, _m, sigma,
			r, key, digest);
	proof_t pf = { y, w, _y, ys, d, p, s, digest };
	serial_writer_init(&bw, published, sizeof(published));
	proof_walk(&pf, &bw, NULL);
	len = bw.overflow ? 0 : serial_bytes(&bw);

	for (int i = 0; i < MSGS; i++) {
		commit_free(&d[i]);
		commit_free(&p[i]);
		nmod_poly_clear(s[i]);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(t[i][k]);
			nmod_poly_clear(tp[i][k]);
			nmod_poly_clear(_t[i][k]);
			nmod_poly_clear(ys[i][0][k]);
			nmod_poly_clear(vs[i][k]);
			nmod_poly_clear(u[i][k]);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_clear(y[i][j][k]);
				nmod_poly_clear(w[i][j][k]);
				nmod_poly_clear(_y[i][j][k]);
			}
		}
	}
	return len;
}

/* Anyone: verify the published proof, from its bytes alone, against the
 * commitments and the published votes. */
static int verify(const uint8_t *buf, size_t len, commit_t com[MSGS],
		nmod_poly_t _m[MSGS], commitkey_t *key) {
	static commit_t d[MSGS], p[MSGS];
	static nmod_poly_t s[MSGS], u[MSGS][2];
	static nmod_poly_t y[MSGS][WIDTH][2], w[MSGS][WIDTH][2], _y[MSGS][WIDTH][2];
	static nmod_poly_t ys[MSGS][1][2], vs[MSGS][2];
	static nmod_poly_t t[MSGS][2], tp[MSGS][2], _t[MSGS][2];
	static uint8_t digest[MSGS][SHA256HashSize];
	bitreader_t br;
	int ok;

	for (int i = 0; i < MSGS; i++) {
		commit_init(&d[i]);
		commit_init(&p[i]);
		nmod_poly_init(s[i], MODP);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_init(t[i][k], MODP);
			nmod_poly_init(tp[i][k], MODP);
			nmod_poly_init(_t[i][k], MODP);
			nmod_poly_init(ys[i][0][k], MODP);
			nmod_poly_init(vs[i][k], MODP);
			nmod_poly_init(u[i][k], MODP);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_init(y[i][j][k], MODP);
				nmod_poly_init(w[i][j][k], MODP);
				nmod_poly_init(_y[i][j][k], MODP);
			}
		}
	}
	memset(digest, 0, sizeof(digest));

	proof_t pf = { y, w, _y, ys, d, p, s, digest };
	serial_reader_init(&br, buf, len);
	proof_walk(&pf, NULL, &br);
	serial_reader_finish(&br);
	ok = !br.overflow && shuffle_verifier(y, w, _y, ys, t, tp, _t, vs, u, d,
			s, com, p, _m, key, digest);

	for (int i = 0; i < MSGS; i++) {
		commit_free(&d[i]);
		commit_free(&p[i]);
		nmod_poly_clear(s[i]);
		for (int k = 0; k < NCRT; k++) {
			nmod_poly_clear(t[i][k]);
			nmod_poly_clear(tp[i][k]);
			nmod_poly_clear(_t[i][k]);
			nmod_poly_clear(ys[i][0][k]);
			nmod_poly_clear(vs[i][k]);
			nmod_poly_clear(u[i][k]);
			for (int j = 0; j < WIDTH; j++) {
				nmod_poly_clear(y[i][j][k]);
				nmod_poly_clear(w[i][j][k]);
				nmod_poly_clear(_y[i][j][k]);
			}
		}
	}
	return ok;
}

static int cmp_int(const void *a, const void *b) {
	return *(const int *) a - *(const int *) b;
}

static void test(flint_rand_t rand) {
	static commitkey_t key;
	static publickey_t pk;
	static privatekey_t sk;
	static commit_t com[MSGS];
	static veritext_t ballot[MSGS];
	static nmod_poly_t vote[MSGS], m[MSGS], _m[MSGS], sigma[MSGS];
	static nmod_poly_t r[MSGS][WIDTH][2];
	static vproof_t vproof[MSGS];
	int pi[MSGS], cast_ok = 1, board_ok = 1, open_ok = 1, same = 1;
	int votes[MSGS], tally[MSGS], count[NCAND] = { 0 };
	size_t len = 0, vbytes = 0;
	uint8_t eid[16], ctx[VCTX_BYTES];
	/* Wall-clock times of the phases, reported at the end; the attack tests
	 * are not counted. */
	double t0, t_key, t_cast, t_board = 0, t_open = 0, t_prove = 0;
	double t_verify = 0;

	t0 = now_s();
	commit_keyinit(&key);
	commit_keygen(&key);
	encrypt_keyinit(&pk, &sk);
	encrypt_keygen(&pk, &sk);
	t_key = now_s() - t0;
	validity_setup();
	randombytes(eid, sizeof(eid));
	for (int i = 0; i < MSGS; i++) {
		vproof[i].proof = (uint8_t *) malloc(validity_proof_len());
		vproof[i].len = 0;
		commit_init(&com[i]);
		vericrypt_init(&ballot[i]);
		nmod_poly_init(vote[i], MODP);
		nmod_poly_init(m[i], MODP);
		nmod_poly_init(_m[i], MODP);
		nmod_poly_init(sigma[i], MODP);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_init(r[i][j][k], MODP);
			}
		}
	}

	/* Voting: the votes are the test oracle, kept apart from what the board,
	 * the authority and the verifier use. */
	t0 = now_s();
	for (int i = 0; i < MSGS; i++) {
		votes[i] = uniform_below(NCAND);
		int_to_monomial(vote[i], votes[i]);
		validity_ctx(ctx, eid, i);
		cast_ok &= cast(&com[i], &ballot[i], &vproof[i], vote[i], &key, &pk, ctx,
				rand);
		vbytes += vproof[i].len;
	}
	t_cast = now_s() - t0;

	TEST_ONCE("every voter casts a ballot") {
		TEST_ASSERT(cast_ok, end);
	} TEST_END;

	TEST_ONCE("the board accepts every ballot against its commitment") {
		t0 = now_s();
		for (int i = 0; i < MSGS; i++) {
			validity_ctx(ctx, eid, i);
			board_ok &= check_ballot(&ballot[i], &vproof[i], &com[i], &key, &pk,
					ctx);
		}
		t_board = now_s() - t0;
		TEST_ASSERT(board_ok, end);
	} TEST_END;

	TEST_ONCE("a ballot copied onto another voter's commitment is rejected") {
		/* A voter cannot reuse another's ballot for its own commitment: the
		 * proof is for the relation of the commitment it was made with. */
		validity_ctx(ctx, eid, 1);
		int copied = check_ballot(&ballot[0], &vproof[0], &com[1], &key, &pk,
				ctx);
		TEST_ASSERT(copied == 0, end);
	} TEST_END;

	TEST_ONCE("a validity proof does not verify for another voter") {
		/* The proof is bound to the voter: moved, with its commitment, to
		 * another voter, it fails. */
		validity_ctx(ctx, eid, 1);
		TEST_ASSERT(check_validity(&vproof[0], &com[0], &key, ctx) == 0, end);
	} TEST_END;

	TEST_ONCE("a vote for two candidates cannot be cast") {
		/* X^0 + X^1, committed honestly: binary, but not one vote. */
		static commit_t c;
		static veritext_t b;
		vproof_t v;
		nmod_poly_t two;

		commit_init(&c);
		vericrypt_init(&b);
		nmod_poly_init(two, MODP);
		nmod_poly_set_coeff_ui(two, 0, 1);
		nmod_poly_set_coeff_ui(two, 1, 1);
		v.proof = (uint8_t *) malloc(validity_proof_len());
		v.len = 0;
		validity_ctx(ctx, eid, 0);
		int ok = cast(&c, &b, &v, two, &key, &pk, ctx, rand);
		free(v.proof);
		nmod_poly_clear(two);
		vericrypt_free(&b);
		commit_free(&c);
		TEST_ASSERT(ok == 0, end);
	} TEST_END;

	TEST_ONCE("the authority opens every commitment from its ballot") {
		t0 = now_s();
		for (int i = 0; i < MSGS; i++) {
			open_ok &= open_ballot(m[i], r[i], &ballot[i], &com[i], &key, &pk,
					&sk);
		}
		t_open = now_s() - t0;
		TEST_ASSERT(open_ok, end);
		/* The oracle: what the authority read is what each voter voted. */
		for (int i = 0; i < MSGS; i++) {
			same &= nmod_poly_equal(m[i], vote[i]);
		}
		TEST_ASSERT(same, end);
	} TEST_END;

	/* The authority shuffles: _m[i] = m[pi(i)], sigma_i = g(pi(i)). */
	for (int i = 0; i < MSGS; i++) {
		pi[i] = i;
	}
	for (int i = MSGS - 1; i > 0; i--) {
		int j = uniform_below(i + 1), x = pi[i];
		pi[i] = pi[j];
		pi[j] = x;
	}
	for (int i = 0; i < MSGS; i++) {
		nmod_poly_set(_m[i], m[pi[i]]);
		int_to_monomial(sigma[i], pi[i]);
	}

	TEST_ONCE("the shuffle proof verifies from its published bytes") {
		t0 = now_s();
		len = prove(com, m, _m, sigma, r, &key);
		t_prove = now_s() - t0;
		TEST_ASSERT(len > 0, end);
		t0 = now_s();
		TEST_ASSERT(verify(published, len, com, _m, &key) == 1, end);
		t_verify = now_s() - t0;
	} TEST_END;

	TEST_ONCE("the tally of the published votes is the votes cast") {
		int valid = 1;

		for (int i = 0; i < MSGS; i++) {
			tally[i] = candidate(_m[i]);
			valid &= tally[i] >= 0;
		}
		TEST_ASSERT(valid, end);
		qsort(tally, MSGS, sizeof(int), cmp_int);
		qsort(votes, MSGS, sizeof(int), cmp_int);
		TEST_ASSERT(memcmp(tally, votes, sizeof(votes)) == 0, end);
	} TEST_END;

	TEST_ONCE("a changed published vote is rejected") {
		/* Count a vote for the next candidate instead, with the same proof. */
		int j = candidate(_m[0]);
		nmod_poly_t keep;

		nmod_poly_init(keep, MODP);
		nmod_poly_set(keep, _m[0]);
		int_to_monomial(_m[0], (j + 1) % NCAND);
		int ok = verify(published, len, com, _m, &key);
		nmod_poly_set(_m[0], keep);
		nmod_poly_clear(keep);
		TEST_ASSERT(ok == 0, end);
	} TEST_END;

	TEST_ONCE("a corrupted proof is rejected") {
		static uint8_t bad[PUBLISHED_CAP];

		memcpy(bad, published, len);
		bad[len / 2] ^= 1;
		TEST_ASSERT(verify(bad, len, com, _m, &key) == 0, end);
	} TEST_END;

	TEST_ONCE("a substituted vote cannot be proven") {
		/* The authority replaces one vote by another candidate and proves the
		 * result with the honest openings and permutation. */
		static nmod_poly_t sub[MSGS];
		size_t n;
		int ok;

		for (int i = 0; i < MSGS; i++) {
			nmod_poly_init(sub[i], MODP);
			nmod_poly_set(sub[i], _m[i]);
		}
		int_to_monomial(sub[0], (candidate(_m[0]) + 1) % NCAND);
		n = prove(com, m, sub, sigma, r, &key);
		ok = n > 0 && verify(published, n, com, sub, &key);
		for (int i = 0; i < MSGS; i++) {
			nmod_poly_clear(sub[i]);
		}
		TEST_ASSERT(ok == 0, end);
	} TEST_END;

	for (int i = 0; i < MSGS; i++) {
		int j = candidate(_m[i]);
		if (j >= 0) {
			count[j]++;
		}
	}
	printf("\n    %d voters, %d candidates; tally:", MSGS, NCAND);
	for (int j = 0; j < NCAND; j++) {
		printf(" %d", count[j]);
	}
	printf("\n    shuffle proof %zu bytes, %.1f KB per vote\n", len,
			len / 1024.0 / MSGS);
	printf("    validity proof %.1f KB on average\n", vbytes / 1024.0 / MSGS);
	printf("    timing, wall clock: key generation %.2f s; casting (and proving "
			"the vote valid) %.2f s, %.2f ms per ballot; board check %.2f s, %.2f "
			"ms per ballot\n",
			t_key, t_cast, 1e3 * t_cast / MSGS, t_board,
			1e3 * t_board / MSGS);
	printf("    timing: authority opens the ballots %.2f s, %.2f ms per "
			"ballot; shuffle prove %.2f s, verify %.2f s\n", t_open,
			1e3 * t_open / MSGS, t_prove, t_verify);

  end:
	for (int i = 0; i < MSGS; i++) {
		free(vproof[i].proof);
		commit_free(&com[i]);
		vericrypt_free(&ballot[i]);
		nmod_poly_clear(vote[i]);
		nmod_poly_clear(m[i]);
		nmod_poly_clear(_m[i]);
		nmod_poly_clear(sigma[i]);
		for (int j = 0; j < WIDTH; j++) {
			for (int k = 0; k < NCRT; k++) {
				nmod_poly_clear(r[i][j][k]);
			}
		}
	}
	encrypt_keyfree(&pk, &sk);
	commit_keyfree(&key);
}

int main(int argc, char *argv[]) {
	flint_rand_t rand;

	setvbuf(stdout, NULL, _IONBF, 0);
	flint_rand_init(rand);
	commit_setup();
	encrypt_setup();

	printf("\n** End-to-end test of the voting protocol:\n\n");
	test(rand);

	encrypt_finish();
	commit_finish();
	flint_rand_clear(rand);
	return test_failures() != 0;
}
