/* Adapted from Intel® Advanced Encryption Standard (Intel® AES) Instructions Set - Rev 3.01
 * https://software.intel.com/sites/default/files/article/165683/aes-wp-2012-09-22-v01.pdf
 */

#include "fastrandombytes.h"
#include "randombytes.h"
#include <string.h>
#include <x86intrin.h>

static prg_t global;
static int global_seeded = 0;
/* The counter is incremented as a 64-bit lane, so it wraps only after 2^64
 * blocks; a 32-bit add would repeat the keystream after 2^32. */
static const __m128i ONE = { 1, 0 };

static inline void KEY_256_ASSIST_1(__m128i *temp1, __m128i *temp2) {
	__m128i temp4;
	*temp2 = _mm_shuffle_epi32(*temp2, 0xff);
	temp4 = _mm_slli_si128(*temp1, 0x4);
	*temp1 = _mm_xor_si128(*temp1, temp4);
	temp4 = _mm_slli_si128(temp4, 0x4);
	*temp1 = _mm_xor_si128(*temp1, temp4);
	temp4 = _mm_slli_si128(temp4, 0x4);
	*temp1 = _mm_xor_si128(*temp1, temp4);
	*temp1 = _mm_xor_si128(*temp1, *temp2);
}

static inline void KEY_256_ASSIST_2(__m128i *temp1, __m128i *temp3) {
	__m128i temp2, temp4;
	temp4 = _mm_aeskeygenassist_si128(*temp1, 0x0);
	temp2 = _mm_shuffle_epi32(temp4, 0xaa);
	temp4 = _mm_slli_si128(*temp3, 0x4);
	*temp3 = _mm_xor_si128(*temp3, temp4);
	temp4 = _mm_slli_si128(temp4, 0x4);
	*temp3 = _mm_xor_si128(*temp3, temp4);
	temp4 = _mm_slli_si128(temp4, 0x4);
	*temp3 = _mm_xor_si128(*temp3, temp4);
	*temp3 = _mm_xor_si128(*temp3, temp2);
}

/* round_key <-- aes256_key_expansion(seed), iv <-- 0 */
void prg_setseed(prg_t *prg, const unsigned char *seed) {
	__m128i *round_key = prg->round_key;
	__m128i temp1, temp2, temp3;

	temp1 = _mm_loadu_si128((__m128i *) seed);
	temp3 = _mm_loadu_si128((__m128i *) (seed + 16));
	round_key[0] = temp1;
	round_key[1] = temp3;
	temp2 = _mm_aeskeygenassist_si128(temp3, 0x01);
	KEY_256_ASSIST_1(&temp1, &temp2);
	round_key[2] = temp1;
	KEY_256_ASSIST_2(&temp1, &temp3);
	round_key[3] = temp3;
	temp2 = _mm_aeskeygenassist_si128(temp3, 0x02);
	KEY_256_ASSIST_1(&temp1, &temp2);
	round_key[4] = temp1;
	KEY_256_ASSIST_2(&temp1, &temp3);
	round_key[5] = temp3;
	temp2 = _mm_aeskeygenassist_si128(temp3, 0x04);
	KEY_256_ASSIST_1(&temp1, &temp2);
	round_key[6] = temp1;
	KEY_256_ASSIST_2(&temp1, &temp3);
	round_key[7] = temp3;
	temp2 = _mm_aeskeygenassist_si128(temp3, 0x08);
	KEY_256_ASSIST_1(&temp1, &temp2);
	round_key[8] = temp1;
	KEY_256_ASSIST_2(&temp1, &temp3);
	round_key[9] = temp3;
	temp2 = _mm_aeskeygenassist_si128(temp3, 0x10);
	KEY_256_ASSIST_1(&temp1, &temp2);
	round_key[10] = temp1;
	KEY_256_ASSIST_2(&temp1, &temp3);
	round_key[11] = temp3;
	temp2 = _mm_aeskeygenassist_si128(temp3, 0x20);
	KEY_256_ASSIST_1(&temp1, &temp2);
	round_key[12] = temp1;
	KEY_256_ASSIST_2(&temp1, &temp3);
	round_key[13] = temp3;
	temp2 = _mm_aeskeygenassist_si128(temp3, 0x40);
	KEY_256_ASSIST_1(&temp1, &temp2);
	round_key[14] = temp1;

	prg->iv = _mm_setzero_si128();
}

static inline void AES_ctr_round(prg_t *prg, unsigned char *out) {
	const __m128i *round_key = prg->round_key;
	__m128i tmp;

	tmp = _mm_xor_si128(prg->iv, round_key[0]);
	tmp = _mm_aesenc_si128(tmp, round_key[1]);
	tmp = _mm_aesenc_si128(tmp, round_key[2]);
	tmp = _mm_aesenc_si128(tmp, round_key[3]);
	tmp = _mm_aesenc_si128(tmp, round_key[4]);
	tmp = _mm_aesenc_si128(tmp, round_key[5]);
	tmp = _mm_aesenc_si128(tmp, round_key[6]);
	tmp = _mm_aesenc_si128(tmp, round_key[7]);
	tmp = _mm_aesenc_si128(tmp, round_key[8]);
	tmp = _mm_aesenc_si128(tmp, round_key[9]);
	tmp = _mm_aesenc_si128(tmp, round_key[10]);
	tmp = _mm_aesenc_si128(tmp, round_key[11]);
	tmp = _mm_aesenc_si128(tmp, round_key[12]);
	tmp = _mm_aesenc_si128(tmp, round_key[13]);
	tmp = _mm_aesenclast_si128(tmp, round_key[14]);
	_mm_storeu_si128((__m128i *) out, tmp);

	prg->iv = _mm_add_epi64(prg->iv, ONE);
}

/* r <-- aes256_ctr(round_key, iv, rlen) */
void prg_bytes(prg_t *prg, unsigned char *r, unsigned long long rlen) {
	unsigned char ct[16];
	unsigned long long num_of_blocks = rlen >> 4;
	unsigned long long i;

	for (i = 0; i < num_of_blocks; i++) {
		AES_ctr_round(prg, r + (i << 4));
	}

	if (rlen & 0x0f) {
		AES_ctr_round(prg, ct);

		memcpy(r + (i << 4), ct, rlen & 0x0f);
	}
}

void fastrandombytes_setseed(const unsigned char *randomness) {
	prg_setseed(&global, randomness);
	global_seeded = 1;
}

void fastrandombytes(unsigned char *r, unsigned long long rlen) {
	if (!global_seeded) {
		unsigned char seed[32];

		randombytes(seed, sizeof(seed));
		fastrandombytes_setseed(seed);
		memset(seed, 0, sizeof(seed));
	}
	prg_bytes(&global, r, rlen);
}
