/* ****************************** *
 * Implemented by Raymond K. ZHAO *
 *                                *
 * PRG                            *
 * ****************************** */
 
#ifndef FASTRANDOMBYTES_H
#define FASTRANDOMBYTES_H

void fastrandombytes(unsigned char *r, unsigned long long rlen);
void fastrandombytes_setseed(const unsigned char *randomness);

#endif
