/* ****************************** *
 * Implemented by Raymond K. ZHAO *
 *                                *
 * Discrete Gaussian Sampler      *
 * ****************************** */
 
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "vcl/vectorclass.h"
#include "vcl/vectormath_exp.h"
#include "vcl/vectormath_trig.h"

#include "randombytes.h"
#include "fastrandombytes.h"
#include "cpucycles.h"
#include "param.h"
#include "gaussian.h"

#include <x86intrin.h>

#define BOX_MULLER_BYTES (8 * 8)

/* Box-Muller sampler precision: 53 bits */
#define FP_PRECISION 53
static const double INT64_DOUBLE_FACTOR = 1.0 / (1LL << FP_PRECISION);
#define FP_MASK ((1LL << FP_PRECISION) - 1)

/* Parameters used by FACCT */
#define COMP_ENTRY_SIZE 9
#define EXP_MANTISSA_PRECISION 52
#define EXP_MANTISSA_MASK ((1LL << EXP_MANTISSA_PRECISION) - 1)
#define R_MANTISSA_PRECISION (EXP_MANTISSA_PRECISION + 1)
#define R_MANTISSA_MASK ((1LL << R_MANTISSA_PRECISION) - 1)
#define R_EXPONENT_L (8 * COMP_ENTRY_SIZE - R_MANTISSA_PRECISION)
#define DOUBLE_ONE (1023LL << 52)

#define NORM_BATCH 8 /* The AVX2 Box-Muller sampler returns 8 samples in batch */
#define DISCRETE_BYTES (1 + NORM_BATCH * COMP_ENTRY_SIZE)

static const double SIGMA = SIGMA_PARAM; /* sigma=2, 32768, 131072 in our benchmark */
static const double DISCRETE_NORMALISATION = 0.5 * M_SQRT1_2 * M_2_SQRTPI * (1.0 / SIGMA); /* 1/S=1/(sigma*sqrt(2*pi)) */
static const double SIGMA_INV = -1.0 / (2 * SIGMA * SIGMA); /* -1/sigma^2 */
static const double M_1_LN2 = 1.0 / M_LN2;

static const __m128d V_SIGMA_INV = {SIGMA_INV, 0};
static const __m128d V_DISCRETE_NORMALISATION = {DISCRETE_NORMALISATION, 0};
static const __m128d V_M_1_LN2 = {M_1_LN2, 0};
static const __m128d V_M_LN2 = {M_LN2, 0};
static const __m128d V_M_1_2 = {0.5, 0};
static const __m128d V_M_1_2_NEG = {-0.5, 0};
static const __m128d V_M_1 = {1.0, 0};
static const __m128d V_M_1_NEG = {-1.0, 0};

#define T 1000
#define TOTAL 1024

static double norm[NORM_BATCH];
static uint64_t head = NORM_BATCH;

/* coefficients of the FACCT exp(x) evaluation polynomial */
static const uint64_t EXP_COFF[] = {0x3e9e44a0d1ea10b4,
									0x3ec36ed25eaa5962,
									0x3efaa0f31065471a,
									0x3f29f143793315db,
									0x3f56c27c690d439b,
									0x3f811105d90e4b7e,
									0x3fa55555e692a339,
									0x3fc555555124dad1,
									0x3fe00000000f3594,
									0x3fefffffffffb3f3,
									0x3ff0000000000000};								   

/* Box-Muller sampler */
static inline void box_muller()
{
	unsigned char r[BOX_MULLER_BYTES];
	Vec4q l1, l2;
	Vec4d r1, r2;
	Vec4d r2_sin, r2_cos;
	
	fastrandombytes(r, BOX_MULLER_BYTES);
	
	l1.load((uint64_t *)r);
	l2.load(((uint64_t *)r) + 4);
	
	r1 = to_double((l1 & FP_MASK) + 1) * INT64_DOUBLE_FACTOR;
	r2 = to_double((l2 & FP_MASK) + 1) * INT64_DOUBLE_FACTOR;
	
	r1 = sqrt(-2.0 * log(r1)) * SIGMA;
	r2 = 2.0 * M_PI * r2;
	
	r2_sin = sincos(&r2_cos, r2);
	
	r2_cos = r1 * r2_cos;
	r2_sin = r1 * r2_sin;
	
	r2_cos.store(norm);
	r2_sin.store(norm + 4);
}

/* FACCT exp(x) */
static inline __m128d exp_constant(const __m128d x)
{
	__m128d vx, vx1, vx2, vsum, vres;
	__m128i vt;
	
	vx = _mm_mul_sd(x, V_SIGMA_INV);

	/* exp(x)=exp(floor(x/ln2)*ln2+a)=2^(floor(x/ln2))*exp(a), where a is in [0,ln2]
	 * we only evaluate exp(a) by using a polynomial */	
	vx1 = _mm_floor_pd(_mm_mul_sd(vx, V_M_1_LN2));
	vt = _mm_slli_epi64(_mm_cvtpd_epi32(vx1), 52);
	vx2 = _mm_sub_sd(vx, _mm_mul_sd(vx1, V_M_LN2));

	/* evaluate 2^a */
	vsum = _mm_add_sd(_mm_mul_sd(_mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[0])), vx2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[1])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[2])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[3])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[4])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[5])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[6])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[7])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[8])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[9])));
	vsum = _mm_add_sd(_mm_mul_sd(vsum, vx2), _mm_castsi128_pd(_mm_cvtsi64x_si128(EXP_COFF[10])));	
	
	/* combine to compute exp(x) */
	vres = _mm_castsi128_pd(_mm_add_epi64(vt, _mm_castpd_si128(vsum)));
	
	return vres;
}

/* FACCT rejection check */
static inline int64_t comp(const unsigned char *r, const uint64_t res)
{
	uint64_t res_mantissa, res_exponent;
	uint64_t r1, r2;
	uint64_t r_mantissa, r_exponent;
	
	res_mantissa = (res & EXP_MANTISSA_MASK) | (1LL << EXP_MANTISSA_PRECISION);
	res_exponent = R_EXPONENT_L - 1023 + 1 + (res >> EXP_MANTISSA_PRECISION); 
	
	r1 = *((uint64_t *)r);
	r2 = (uint64_t)(r[8]);
	
	r_mantissa = r1 & R_MANTISSA_MASK;
	r_exponent = (r1 >> R_MANTISSA_PRECISION) | (r2 << (64 - R_MANTISSA_PRECISION));
	
	return (((1LL << 63) ^ ((res - DOUBLE_ONE) | (DOUBLE_ONE - res))) | ((r_mantissa - res_mantissa) & (r_exponent - (1LL << res_exponent)))) >> 63;	
}

/* We can merge branches in Algorithm 2 in the paper as follows: 
 * Let y<--N(0, sigma^2)
 * Let b<--U({0,1})
 * If (b==0):
 *     Let y_r=round(y)-1
 *     Let boolean cmp=(y<=0.5)
 * Else:
 *     Let y_r=round(y)+1
 *     Let boolean cmp=(y>=-0.5)
 * EndIf
 * If (cmp is true):
 *     Let r<--U([0,1))
 *     If (r<exp(-((y_r+c_F)^2-y^2)/2sigma^2):
 *         Return y_r
 *     Else:
 *         Restart
 *     EndIf
 * Else:
 *     Restart
 * EndIf
 */
int64_t discrete_gaussian(const double center)
{
	unsigned char r[DISCRETE_BYTES]; 	
	
	uint64_t b, i, dcmp, dcmp1, drej, drc, dcmprc;
	
	__m128d vc, vcr, vy, vyr, vyrc;
	__m128i vb, vcmprc; 
	
	vc = _mm_load_sd(&center);
	vcr = _mm_round_pd(vc, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC);
	vc = _mm_sub_sd(vcr, vc);
	
	_mm_storel_epi64((__m128i *)(&drc), _mm_castpd_si128(_mm_mul_sd(exp_constant(_mm_mul_sd(vc, vc)), V_DISCRETE_NORMALISATION)));

	fastrandombytes(r, COMP_ENTRY_SIZE);
	dcmprc = -comp(r, drc);
	vcmprc = _mm_loadl_epi64((__m128i *)(&dcmprc));

	while (true)
	{
		fastrandombytes(r, DISCRETE_BYTES);
		
		for (i = 0; i < 8; i++, head++)
		{
			if (head >= NORM_BATCH)
			{
				head = 0;
				box_muller();
			}

			vy = _mm_load_sd(norm + head);
			b = (r[DISCRETE_BYTES - 1] >> i) & 0x01;
						
			b = -((-b | b) >> 63);
			dcmp = _mm_comile_sd(vy, V_M_1_2);
			dcmp1 = dcmp ^ _mm_comige_sd(vy, V_M_1_2_NEG);
			dcmp ^= b & dcmp1;

			vb = _mm_loadl_epi64((__m128i *)(&b));
			vyr = _mm_add_sd(_mm_round_pd(vy, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC), _mm_castsi128_pd(_mm_xor_si128(_mm_castpd_si128(V_M_1_NEG), _mm_and_si128(vb, _mm_xor_si128(_mm_castpd_si128(V_M_1_NEG), _mm_castpd_si128(V_M_1))))));

			if (dcmp)
			{
				vyrc = _mm_add_sd(vyr, vc);
				_mm_storel_epi64((__m128i *)(&drej), _mm_castpd_si128(exp_constant(_mm_mul_sd(_mm_add_sd(vyrc, vy), _mm_sub_sd(vyrc, vy)))));

				if (comp(r + i * COMP_ENTRY_SIZE, drej))
				{
					head++;
					return _mm_cvtsd_si64(_mm_add_sd(_mm_castsi128_pd(_mm_xor_si128(_mm_castpd_si128(vyr), _mm_and_si128(vcmprc, _mm_castpd_si128(vyr)))), vcr));
				}				
			}
		}
	}
}

#ifdef MAIN
int main()
{
	unsigned char seed[32];
	uint64_t i, t;
	int64_t sample[TOTAL];
	double center[TOTAL];
	
	long long cycle1, cycle2;
	srand(time(NULL));	
for (t = 0; t < T; t++)
{
	randombytes(seed, 32);
	
	fastrandombytes_setseed(seed);
	
	head = NORM_BATCH;

	for (i = 0; i < TOTAL; i++)
	{
		center[i] = (double)(rand()) / (double)(RAND_MAX);
	}
	
	cycle1 = cpucycles();
	for (i = 0; i < TOTAL; i++)
	{
		sample[i] = discrete_gaussian(center[i]);
	}
	cycle2 = cpucycles();
	
	for (i = 0; i < TOTAL; i++)
	{
		printf("%lld ", sample[i]);
	}
	printf("\n");
	printf("cycle:%lld\n", cycle2 - cycle1);
}
}
#endif
