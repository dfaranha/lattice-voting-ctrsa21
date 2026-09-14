/**
 * @file
 *
 * Implementation of the LNP-style proof machinery for the is_bin sub-proof.
 *
 * @ingroup lnp
 */

#include <assert.h>

#include "param.h"
#include "lnp.h"
#include "test.h"
#include "bench.h"
#include "fastrandombytes.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/* The automorphism sigma_k maps X^i to X^(k*i). Exponents live modulo 2*DEGREE
 * because X^DEGREE = -1, so an exponent that lands in [DEGREE, 2*DEGREE) comes
 * back with a sign flip. */
static void auto_exponent(slong *e, int *neg, slong i, slong k) {
	slong t = (i % (2 * DEGREE)) * (k % (2 * DEGREE)) % (2 * DEGREE);

	if (t < 0) {
		t += 2 * DEGREE;
	}
	*neg = (t >= DEGREE);
	*e = *neg ? t - DEGREE : t;
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

int lnp_auto_swaps(slong k) {
	slong r = k % 4;

	if (r < 0) {
		r += 4;
	}
	/* X^(DEGREE/2) squares to -1, so raising it to an odd power k gives
	 * X^(DEGREE/2) when k = 1 mod 4 and -X^(DEGREE/2) when k = 3 mod 4. Since
	 * the two factors are X^(DEGREE/2) + P0 and X^(DEGREE/2) - P0, the second
	 * case sends one factor onto the other. */
	assert(r == 1 || r == 3);
	return (r == 3);
}

void lnp_auto(nmod_poly_t c, const nmod_poly_t a, slong k) {
	nmod_poly_t t;
	slong e;
	int neg;

	assert(k % 2 != 0);

	nmod_poly_init(t, MODP);
	nmod_poly_fit_length(t, DEGREE);
	for (slong i = 0; i < DEGREE; i++) {
		ulong v = nmod_poly_get_coeff_ui(a, i);

		if (v == 0) {
			continue;
		}
		auto_exponent(&e, &neg, i, k);
		if (neg) {
			v = nmod_neg(v, a->mod);
		}
		/* sigma_k is a bijection on exponents, so no two terms collide and a
		 * plain store would do; adding keeps this correct if a ever carries
		 * more than DEGREE coefficients. */
		nmod_poly_set_coeff_ui(t, e,
				nmod_add(nmod_poly_get_coeff_ui(t, e), v, a->mod));
	}
	nmod_poly_set(c, t);
	nmod_poly_clear(t);
}

void lnp_auto_crt(pcrt_poly_t c, pcrt_poly_t a, slong k) {
	nmod_poly_t t[NCRT];
	int swap = lnp_auto_swaps(k);

	/* The residue of a modulo one factor determines the residue of sigma_k(a)
	 * modulo the image factor, because sigma_k maps the one ideal onto the
	 * other. So it suffices to apply sigma_k to the canonical lift of each
	 * component and reduce the image into the target component. */
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(t[i], MODP);
		lnp_auto(t[i], a[i], k);
	}
	for (int i = 0; i < NCRT; i++) {
		int j = swap ? (NCRT - 1 - i) : i;

		pcrt_poly_reduce(c[j], t[i], j);
	}
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(t[i]);
	}
}

/*============================================================================*/
/* Tests                                                                      */
/*============================================================================*/

/* sigma_{-1} in coefficient representation, written out directly from the
 * definition as an independent check on the general routine. */
static void auto_inverse_ref(nmod_poly_t c, nmod_poly_t a) {
	nmod_poly_t t;

	nmod_poly_init(t, MODP);
	nmod_poly_fit_length(t, DEGREE);
	nmod_poly_set_coeff_ui(t, 0, nmod_poly_get_coeff_ui(a, 0));
	for (slong i = 1; i < DEGREE; i++) {
		/* X^-i = -X^(DEGREE-i) because X^DEGREE = -1. */
		nmod_poly_set_coeff_ui(t, DEGREE - i,
				nmod_neg(nmod_poly_get_coeff_ui(a, i), a->mod));
	}
	nmod_poly_set(c, t);
	nmod_poly_clear(t);
}

static void test_auto(flint_rand_t rng) {
	nmod_poly_t a, b, d, e;
	pcrt_poly_t ca, cb;
	/* 2*DEGREE-1 = -1 mod 2*DEGREE, so this is sigma_{-1}. */
	const slong minus1 = 2 * DEGREE - 1;

	nmod_poly_init(a, MODP);
	nmod_poly_init(b, MODP);
	nmod_poly_init(d, MODP);
	nmod_poly_init(e, MODP);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_init(ca[i], MODP);
		nmod_poly_init(cb[i], MODP);
	}

	TEST_BEGIN("sigma_{-1} matches its definition") {
		commit_sample_rand(a, rng, DEGREE);
		lnp_auto(b, a, minus1);
		auto_inverse_ref(d, a);
		TEST_ASSERT(nmod_poly_equal(b, d) == 1, end);
	} TEST_END;

	TEST_BEGIN("sigma_k is a ring homomorphism") {
		commit_sample_rand(a, rng, DEGREE);
		commit_sample_rand(b, rng, DEGREE);
		/* sigma_k(a*b) = sigma_k(a)*sigma_k(b). */
		commit_poly_mulmod(d, a, b);
		lnp_auto(d, d, minus1);
		lnp_auto(e, a, minus1);
		lnp_auto(b, b, minus1);
		commit_poly_mulmod(e, e, b);
		TEST_ASSERT(nmod_poly_equal(d, e) == 1, end);
	} TEST_END;

	TEST_BEGIN("sigma_{-1} is an involution") {
		commit_sample_rand(a, rng, DEGREE);
		lnp_auto(b, a, minus1);
		lnp_auto(b, b, minus1);
		TEST_ASSERT(nmod_poly_equal(a, b) == 1, end);
	} TEST_END;

	TEST_BEGIN("constant coefficient of sigma_{-1}(a)*b is <a, b>") {
		commit_sample_rand(a, rng, DEGREE);
		commit_sample_rand(b, rng, DEGREE);
		lnp_auto(d, a, minus1);
		commit_poly_mulmod(d, d, b);
		nmod_poly_zero(e);
		for (slong i = 0; i < DEGREE; i++) {
			nmod_poly_set_coeff_ui(e, 0,
					nmod_add(nmod_poly_get_coeff_ui(e, 0),
					nmod_mul(nmod_poly_get_coeff_ui(a, i),
					nmod_poly_get_coeff_ui(b, i), a->mod), a->mod));
		}
		TEST_ASSERT(nmod_poly_get_coeff_ui(d, 0) ==
				nmod_poly_get_coeff_ui(e, 0), end);
	} TEST_END;

	TEST_BEGIN("CRT-aware sigma_k agrees with the coefficient reference") {
		for (slong k = 1; k < 8; k += 2) {
			commit_sample_rand(a, rng, DEGREE);
			/* Reference: map in coefficient representation, then reduce. */
			lnp_auto(b, a, k);
			pcrt_poly_reduce(cb[0], b, 0);
			pcrt_poly_reduce(cb[1], b, 1);
			/* Under test: reduce first, then map in CRT representation. */
			pcrt_poly_reduce(ca[0], a, 0);
			pcrt_poly_reduce(ca[1], a, 1);
			lnp_auto_crt(ca, ca, k);
			TEST_ASSERT(nmod_poly_equal(ca[0], cb[0]) == 1, end);
			TEST_ASSERT(nmod_poly_equal(ca[1], cb[1]) == 1, end);
		}
	} TEST_END;

	TEST_ONCE("sigma_{-1} exchanges the two CRT factors") {
		/* Applying sigma_{-1} to the first irreducible factor must give a
		 * polynomial that vanishes modulo the second, not modulo the first. */
		lnp_auto(a, *commit_irred(0), minus1);
		pcrt_poly_reduce(ca[0], a, 0);
		pcrt_poly_reduce(ca[1], a, 1);
		TEST_ASSERT(nmod_poly_is_zero(ca[1]) == 1, end);
		TEST_ASSERT(nmod_poly_is_zero(ca[0]) == 0, end);
		TEST_ASSERT(lnp_auto_swaps(minus1) == 1, end);
		TEST_ASSERT(lnp_auto_swaps(5) == 0, end);
	} TEST_END;

  end:
	nmod_poly_clear(a);
	nmod_poly_clear(b);
	nmod_poly_clear(d);
	nmod_poly_clear(e);
	for (int i = 0; i < NCRT; i++) {
		nmod_poly_clear(ca[i]);
		nmod_poly_clear(cb[i]);
	}
}

/* Select which phases to run: "test", "bench", or neither for both. Same
 * helper as in commit.c and shuffle.c; it stays local to each binary so that
 * the test harness does not have to be shared between them. */
static int phase_selected(int argc, char *argv[], const char *phase) {
	return argc < 2 || strcmp(argv[1], phase) == 0;
}

int main(int argc, char *argv[]) {
	flint_rand_t rand;

	flint_rand_init(rand);
	commit_setup();

	if (phase_selected(argc, argv, "test")) {
		printf("\n** Tests for the LNP proof machinery:\n\n");
		test_auto(rand);
	}

	commit_finish();
	flint_rand_clear(rand);
	return 0;
}
