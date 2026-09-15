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
 * That assumption is the number of bits a coefficient needs. A uniform ring
 * element is written at ceil(log2 p) bits per coefficient, and a Gaussian one
 * at ceil(log2 12 sigma), which covers six standard deviations either side.
 * The verifier's norm bound is looser than that: it bounds the whole vector in
 * l2, so a single coefficient could in principle be much larger and still
 * pass. Packing therefore refuses rather than truncating when a coefficient
 * does not fit, which turns the modelling assumption into something an honest
 * run can falsify.
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
 * Bits needed per coefficient.
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
