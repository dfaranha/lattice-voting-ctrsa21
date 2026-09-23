/**
 * @defgroup serial Proof serialisation.
 */
/**
 * @file
 *
 * Bit-level packing for the values a proof transmits, so that its size can be
 * measured rather than modelled. Every size figure in LNP-PARAMS.md used to be
 * computed by hand from the shapes of the structures; this makes the same
 * count come out of the code, and makes the assumption behind it testable.
 *
 * A uniform ring element is written at ceil(log2 p) bits per coefficient,
 * where there is nothing to exploit. A Gaussian one is not: writing it at the
 * width of its bound, ceil(log2 12 sigma), spends the difference between that
 * bound and the distribution's entropy, which is log2(sigma sqrt(2 pi e)). At
 * SIGMA_C that is 20 bits against 17.8. So a Gaussian coefficient is zigzagged
 * and Golomb-Rice coded about its own width instead, which costs about half a
 * bit above the entropy.
 *
 * This also removes an assumption the flat encoding had to make. The verifier
 * bounds the whole vector in l2, so a single coefficient can exceed six
 * standard deviations and still pass; the flat packing had to refuse such a
 * coefficient rather than truncate it. Rice has no width to overflow, and the
 * long-run escape encodes any value at all, so the question does not arise.
 *
 * @ingroup serial
 */

#ifndef SERIAL_H
#define SERIAL_H

#include <stddef.h>
#include <stdint.h>

#include <flint/nmod_poly.h>

#include "param.h"

/*============================================================================*/
/* Type definitions                                                           */
/*============================================================================*/

typedef struct _bitwriter_t {
	uint8_t *buf;
	size_t cap;					/* capacity in bytes */
	size_t bits;				/* bits written so far */
	int overflow;				/* set if a write did not fit */
} bitwriter_t;

typedef struct _bitreader_t {
	const uint8_t *buf;
	size_t cap;
	size_t bits;
	int overflow;
} bitreader_t;

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

void serial_writer_init(bitwriter_t *w, uint8_t *buf, size_t cap);
void serial_reader_init(bitreader_t *r, const uint8_t *buf, size_t cap);

/**
 * Bits needed per coefficient. The Gaussian one is the width a flat encoding
 * would take, which nothing writes at any more; it is what the saving of the
 * entropy coding is measured against.
 */
int serial_bits_uniform(void);
int serial_bits_gauss(ulong sigma);

/**
 * Number of whole bytes written so far, rounding up.
 */
size_t serial_bytes(const bitwriter_t *w);

/**
 * Pack and unpack one ring element held in CRT representation, whose
 * coefficients are uniform. The two components are written as they stand,
 * which is DEGREE coefficients in total.
 */
void serial_put_uniform(bitwriter_t *w, nmod_poly_t a[2]);
void serial_get_uniform(bitreader_t *r, nmod_poly_t a[2]);

/**
 * Pack and unpack one ring element held in CRT representation whose
 * *reconstructed* coefficients are short. The CRT components of a short
 * element are not themselves short, so this reconstructs first and writes the
 * centred coefficients.
 */
void serial_put_gauss(bitwriter_t *w, nmod_poly_t a[2], ulong sigma);
void serial_get_gauss(bitreader_t *r, nmod_poly_t a[2], ulong sigma);

/**
 * Pack and unpack one raw byte, for the Fiat-Shamir digest.
 */
void serial_put_byte(bitwriter_t *w, uint8_t v);
uint8_t serial_get_byte(bitreader_t *r);

#endif /* !SERIAL_H */
