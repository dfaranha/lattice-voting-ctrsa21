/* ****************************** *
 * Implemented by Raymond K. ZHAO *
 *                                *
 * PRG                            *
 * ****************************** */

#ifndef FASTRANDOMBYTES_H
#define FASTRANDOMBYTES_H

#include <x86intrin.h>

/* State of an AES-256-CTR stream. */
typedef struct {
	__m128i round_key[15];
	__m128i iv;
} prg_t;

/* Key a stream from a 32-byte seed and reset its counter. */
void prg_setseed(prg_t *prg, const unsigned char *seed);
/* Read rlen bytes from a stream. */
void prg_bytes(prg_t *prg, unsigned char *r, unsigned long long rlen);

/* The global stream that feeds the Gaussian sampler. It holds the secret masks,
 * so it seeds itself from the operating system on first use and must never be
 * keyed with public data: anything derived from a hash, such as a challenge,
 * belongs in a local prg_t instead. */
void fastrandombytes(unsigned char *r, unsigned long long rlen);
/* Replace the global seed, for reproducible sampler tests only. */
void fastrandombytes_setseed(const unsigned char *randomness);

#endif
