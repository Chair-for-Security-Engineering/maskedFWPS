/* Copyright 2022 UCLouvain, Belgium and PQM4 contributors
 *
 * This file is part of pqm4_masked.
 *
 * pqm4_masked is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, version 3.
 *
 * pqm4_masked is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * pqm4_masked. If not, see <https://www.gnu.org/licenses/>.
 * 
 * Edited by SecEng, RUB, Germany
 */
#include "gadgets.h"
#include "masked.h"
#include "masked_representations.h"
#include "params.h"
#include <stdint.h>

static inline uint32_t pini_and_core(uint32_t a, uint32_t b, uint32_t r) {
#if 0 // original ARM code
  asm("eor %[temp], %[b], %[r]\n\t"
      "and %[temp], %[a], %[temp]\n\t"
      "bic %[s], %[r], %[a]\n\t"
      "eor %[s], %[s], %[temp]"
      : [ s ] "=r"(s),
        [ temp ] "=&r"(
            temp) /* outputs, use temp as an arbitrary-location clobber */
      : [ a ] "r"(a), [ b ] "r"(b), [ r ] "r"(r) /* inputs */
  );
  return s;
#else // C placeholder code
  uint32_t temp = b ^ r;
  temp = temp & a;
  uint32_t s = r&(~a);
  return s^temp;
#endif
}

/*************************************************
 * Name:        masked_and
 *
 * Description: Performs masked AND (z = a & b ) gate with nshares.
 *
 * Arguments:   - size_t nshares: number of shares
 *            - uint32_t *z: output buffer
 *            - size_t z_stride: output buffer stride
 *            - uint32_t *a: first input buffer
 *            - size_t a_stride: a buffer stride
 *            - uint32_t *b: second input buffer
 *            - size_t b_stride: b buffer stride
 **************************************************/
void masked_and(size_t nshares, uint32_t *z, size_t z_stride, const uint32_t *a,
                size_t a_stride, const uint32_t *b, size_t b_stride) {
  uint32_t ztmp[nshares];
  uint32_t r;
  uint32_t i, j;

  for (i = 0; i < nshares; i++) {
    ztmp[i] = a[i * a_stride] & b[i * b_stride];
  }

  for (i = 0; i < (nshares - 1); i++) {
    for (j = i + 1; j < nshares; j++) {
      r = get_random();
      // PINI
      ztmp[i] ^= pini_and_core(a[i * a_stride], b[j * b_stride], r);
      ztmp[j] ^= pini_and_core(a[j * a_stride], b[i * b_stride], r);
    }
  }
  for (i = 0; i < nshares; i++) {
    z[i * z_stride] = ztmp[i];
  }
}


/*************************************************
 * Name:        masked_and_1bit
 *
 * Description: Performs masked AND (z = a & b ) gate with nshares only for the LSB
 *
 * Arguments:   - size_t nshares: number of shares
 *            - uint32_t *z: output buffer
 *            - size_t z_stride: output buffer stride
 *            - uint32_t *a: first input buffer
 *            - size_t a_stride: a buffer stride
 *            - uint32_t *b: second input buffer
 *            - size_t b_stride: b buffer stride
 **************************************************/
void masked_and_1bit(size_t nshares, uint32_t *z, size_t z_stride, const uint32_t *a,
                size_t a_stride, const uint32_t *b, size_t b_stride) {
  uint32_t ztmp[nshares];
  uint32_t i, j;
  static uint8_t rcnt = 0;
  static uint32_t r = 0;

  for (i = 0; i < nshares; i++) {
    ztmp[i] = a[i * a_stride] & b[i * b_stride];
  }

  for (i = 0; i < (nshares - 1); i++) {
    for (j = i + 1; j < nshares; j++) {
      if (rcnt == 0)
      {
        r = get_random();
        rcnt = 32;
      } else {
        r >>= 1;
        rcnt -= 1;
      }
      // PINI
      ztmp[i] ^= pini_and_core(a[i * a_stride], b[j * b_stride], r&1);
      ztmp[j] ^= pini_and_core(a[j * a_stride], b[i * b_stride], r&1);
    }
  }
  for (i = 0; i < nshares; i++) {
    z[i * z_stride] = ztmp[i];
  }
}

/*************************************************
 * Name:        masked_or
 *
 * Description: Performs masked AND (z = a | b ) gate with nshares.
 *              a | b = ~( ~a & ~b )
 *
 * Arguments:   - size_t nshares: number of shares
 *            - uint32_t *z: output buffer
 *            - size_t z_stride: output buffer stride
 *            - uint32_t *a: first input buffer
 *            - size_t a_stride: a buffer stride
 *            - uint32_t *b: second input buffer
 *            - size_t b_stride: b buffer stride
 **************************************************/
void masked_or(size_t nshares, uint32_t *z, size_t z_stride, const uint32_t *a,
                size_t a_stride, const uint32_t *b, size_t b_stride) {
  uint32_t nota[nshares];
  uint32_t notb[nshares];
  uint32_t i;

  nota[0] = ~a[0];
  notb[0] = ~b[0];
  for(i = 1; i<nshares; i++) {
    nota[i] = a[i*a_stride];
    notb[i] = b[i*b_stride];
  }
  
  masked_and(nshares, z, z_stride, nota, 1, notb, 1);
  z[0] = ~z[0];              
}

/*************************************************
 * Name:        masked_or_1bit
 *
 * Description: Performs masked AND (z = a | b ) gate with nshares.
 *              a | b = ~( ~a & ~b )
 *
 * Arguments:   - size_t nshares: number of shares
 *            - uint32_t *z: output buffer
 *            - size_t z_stride: output buffer stride
 *            - uint32_t *a: first input buffer
 *            - size_t a_stride: a buffer stride
 *            - uint32_t *b: second input buffer
 *            - size_t b_stride: b buffer stride
 **************************************************/
void masked_or_1bit(size_t nshares, uint32_t *z, size_t z_stride, const uint32_t *a,
                size_t a_stride, const uint32_t *b, size_t b_stride) {
  uint32_t nota[nshares];
  uint32_t notb[nshares];
  uint32_t i;

  nota[0] = ~a[0];
  notb[0] = ~b[0];
  for(i = 1; i<nshares; i++) {
    nota[i] = a[i*a_stride];
    notb[i] = b[i*b_stride];
  }
  
  masked_and_1bit(nshares, z, z_stride, nota, 1, notb, 1);
  z[0] = ~z[0];              
}

/*************************************************
 * Name:        masked_xor
 *
 * Description: Performs masked XOR (z = a ^ b ) gate with nshares.
 *
 * Arguments:   - size_t nshares: number of shares
 *            - uint32_t *z: output buffer
 *            - size_t z_stride: output buffer stride
 *            - uint32_t *a: first input buffer
 *            - size_t a_stride: a buffer stride
 *            - uint32_t *b: second input buffer
 *            - size_t b_stride: b buffer stride
 **************************************************/
void masked_xor(size_t nshares, uint32_t *out, size_t out_stride,
                const uint32_t *ina, size_t ina_stride, const uint32_t *inb,
                size_t inb_stride) {
  for (size_t i = 0; i < nshares; i++) {
    out[i * out_stride] = ina[i * ina_stride] ^ inb[i * inb_stride];
  }
}

/*************************************************
 * Name:        unmask_boolean
 *
 * Description: Unmask sharing by XORing all the shares
 *
 * Arguments: - size_t nshares: number of shares
 *            - uint32_t *in: input buffer
 *            - size_t in_stride: in buffer stride
 * Returns:   - uint32_t out: masked value
 **************************************************/
uint32_t unmask_boolean(size_t nshares, const uint32_t *in, size_t in_stride) {
  uint32_t out, d;
  out = 0;
  for (d = 0; d < nshares; d++) {
    out ^= in[d * in_stride];
  }
  return out;
}

/*************************************************
 * Name:        copy_sharing
 *
 * Description: Copy input sharing to output sharing 
 *
 * Arguments: - size_t nshares: number of shares
 *            - uint32_t *out: output buffer
 *            - size_t out_stride: out buffer stride
 *            - uint32_t *in: input buffer
 *            - size_t in_stride: in buffer stride
 **************************************************/
void copy_sharing(size_t nshares, uint32_t *out, size_t out_stride,
                  const uint32_t *in, size_t in_stride) {
  for (size_t i = 0; i < nshares; i++) {
    out[i * out_stride] = in[i * in_stride];
  }
}

/*************************************************
 * Name:        secadd
 *
 * Description: Performs masked addition of two bitslice words. 
 *              out = (in1 + in2)%(2**kbits_out) 
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words
 *            - size_t kbits_out: number of bits in the output word.
 *                kbits = kbits_out or kbits = kbits_out - 1
 *            - uint32_t *out: output buffer
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - uint32_t *in1: first input buffer
 *            - size_t in1_msk_stride: stride between shares
 *            - size_t in1_data_stride: stride between masked bits
 *            - uint32_t *in2: second input buffer
 *            - size_t in2_msk_stride: stride between shares
 *            - size_t in2_data_stride: stride between masked bits
 **************************************************/
void secadd(size_t nshares, size_t kbits, size_t kbits_out, uint32_t *out,
            size_t out_msk_stride, size_t out_data_stride, const uint32_t *in1,
            size_t in1_msk_stride, size_t in1_data_stride, const uint32_t *in2,
            size_t in2_msk_stride, size_t in2_data_stride) {
  
  size_t i, d;
  uint32_t carry[nshares];
  uint32_t xpy[nshares];
  uint32_t xpc[nshares];

  masked_and(nshares, carry, 1, &in1[0 * in1_data_stride], in1_msk_stride,
             &in2[0 * in2_data_stride], in2_msk_stride);

  masked_xor(nshares, &out[0 * out_data_stride], out_msk_stride,
             &in1[0 * in1_data_stride], in1_msk_stride,
             &in2[0 * in2_data_stride], in2_msk_stride);

  for (i = 1; i < kbits; i++) {
    // xpy = in2 ^ in1
    // xpc = in1 ^ carry
    // out = xpy ^ carry
    for (d = 0; d < nshares; d++) {
      xpy[d] = in1[i * in1_data_stride + d * in1_msk_stride] ^
               in2[i * in2_data_stride + d * in2_msk_stride];
      xpc[d] = in1[i * in1_data_stride + d * in1_msk_stride] ^ carry[d];
      out[i * out_data_stride + d * out_msk_stride] = xpy[d] ^ carry[d];
    }

    if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
      break;
    } else if (i == (kbits - 1)) {
      masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
      masked_xor(nshares, &out[(kbits)*out_data_stride], out_msk_stride, carry,
                 1, &in1[i * in1_data_stride], in1_msk_stride);
      break;
    }

    masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
    masked_xor(nshares, carry, 1, carry, 1, &in1[i * in1_data_stride],
               in1_msk_stride);
  }

}

/*************************************************
 * Name:        secadd_constant_bmsk
 *
 * Description: Performs masked addition of an input with a masked
 *              constant sucht that:
 *              out = (in + bmsk*constant)%(2**kbits_out) 
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words
 *            - size_t kbits_out: number of bits in the output word.
 *                kbits = kbits_out or kbits = kbits_out - 1
 *            - uint32_t *out: output buffer
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - uint32_t *in: first input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 *            - uint32_t constant: public constant
 *            - uint32_t *bmsk: masked bit buffer
 *            - size_t bmsk_msk_stride: shares stride of masked bit
 **************************************************/
void secadd_constant_bmsk(size_t nshares, size_t kbits, size_t kbits_out,
                          uint32_t *out, size_t out_msk_stride,
                          size_t out_data_stride, const uint32_t *in1,
                          size_t in1_msk_stride, size_t in1_data_stride,
                          uint32_t constant, const uint32_t *bmsk,
                          size_t bmsk_msk_stride) {
  size_t i, d;
  uint32_t carry[nshares];
  uint32_t xpy[nshares];
  uint32_t xpc[nshares];

  if (constant & 0x1) {
    masked_and(nshares, carry, 1, &in1[0 * in1_data_stride], in1_msk_stride,
               bmsk, bmsk_msk_stride);
    masked_xor(nshares, &out[0 * out_data_stride], out_msk_stride,
               &in1[0 * in1_data_stride], in1_msk_stride, bmsk,
               bmsk_msk_stride);
  } else {
    for (d = 0; d < nshares; d++) {
      carry[d] = 0;
    }
    copy_sharing(nshares, out, out_msk_stride, in1, in1_msk_stride);
  }
  for (i = 1; i < kbits; i++) {
    // xpy = in2 ^ in1
    // xpc = in1 ^ carry
    // out = xpy ^ carry
    if ((constant >> i) & 0x1) {
      for (d = 0; d < nshares; d++) {
        xpy[d] = in1[i * in1_data_stride + d * in1_msk_stride] ^
                 bmsk[d * bmsk_msk_stride];
        xpc[d] = in1[i * in1_data_stride + d * in1_msk_stride] ^ carry[d];
        out[i * out_data_stride + d * out_msk_stride] = xpy[d] ^ carry[d];
      }

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
        masked_xor(nshares, &out[(kbits)*out_data_stride], out_msk_stride,
                   carry, 1, &in1[i * in1_data_stride], in1_msk_stride);
        return;
      }

      masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
      masked_xor(nshares, carry, 1, carry, 1, &in1[i * in1_data_stride],
                 in1_msk_stride);
    } else {
      // compute the carry
      masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, carry, 1,
                 &in1[i * in1_data_stride], 1);

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, &out[(kbits)*out_data_stride], out_msk_stride,
                   carry, 1, &in1[i * in1_data_stride], in1_msk_stride);
        return;
      }
      masked_and(nshares, carry, 1, carry, 1, &in1[i * in1_data_stride],
                 in1_msk_stride);
    }
  }
}

/*************************************************
 * Name:        secadd_constant
 *
 * Description: Performs masked addition of an input with a 
 *              public constant such that
 *              out = (in + constant)%(2**kbits_out) 
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words
 *            - size_t kbits_out: number of bits in the output word.
 *                kbits = kbits_out or kbits = kbits_out - 1
 *            - uint32_t *out: output buffer
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - uint32_t *in: first input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 *            - uint32_t constant: public constant
 **************************************************/
void secadd_constant(size_t nshares, size_t kbits, size_t kbits_out,
                     uint32_t *out, size_t out_msk_stride,
                     size_t out_data_stride, const uint32_t *in1,
                     size_t in1_msk_stride, size_t in1_data_stride,
                     uint32_t constant) {

  size_t i, d;
  uint32_t carry[nshares];
  uint32_t xpy[nshares];
  uint32_t xpc[nshares];

  if (constant & 0x1) {
    for (d = 0; d < nshares; d++) {
      carry[d] = in1[d * in1_msk_stride];
    }
    copy_sharing(nshares, out, out_msk_stride, in1, in1_msk_stride);
    out[0] ^= 0xFFFFFFFF;
  } else {
    for (d = 0; d < nshares; d++) {
      carry[d] = 0;
    }
    copy_sharing(nshares, out, out_msk_stride, in1, in1_msk_stride);
  }

  for (i = 1; i < kbits; i++) {
    if ((constant >> i) & 0x1) {
      for (d = 0; d < nshares; d++) {
        xpy[d] = in1[i * in1_data_stride + d * in1_msk_stride];
        xpc[d] = in1[i * in1_data_stride + d * in1_msk_stride] ^ carry[d];
        out[i * out_data_stride + d * out_msk_stride] = xpy[d] ^ carry[d];
      }
      xpy[0] ^= 0xFFFFFFFF;
      out[i * out_data_stride] ^= 0xFFFFFFFF;

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
        masked_xor(nshares, &out[(kbits)*out_data_stride], out_msk_stride,
                   carry, 1, &in1[i * in1_data_stride], in1_msk_stride);

        // add the kbits_out of the constant
        out[(kbits)*out_data_stride] ^=
            0xFFFFFFFF * ((constant >> kbits) & 0x1);

        return;
      }

      masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
      masked_xor(nshares, carry, 1, carry, 1, &in1[i * in1_data_stride],
                 in1_msk_stride);
    } else {
      // compute the carry
      masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, carry, 1,
                 &in1[i * in1_data_stride], 1);

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, &out[(kbits)*out_data_stride], out_msk_stride,
                   carry, 1, &in1[i * in1_data_stride], in1_msk_stride);

        // add the kbits_out of the constant
        out[(kbits)*out_data_stride] ^=
            0xFFFFFFFF * ((constant >> kbits) & 0x1);
        return;
      }
      masked_and(nshares, carry, 1, carry, 1, &in1[i * in1_data_stride],
                 in1_msk_stride);
    }
  }
}

/*************************************************
 * Name:        secadd_modp
 *
 * Description: Performs masked addition of two bitslice words with 
 *              arbitrary modulo such that
 *              out = (in1 + in2)%(p) 
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words:
 *                requires kbits = ceil(log(p))
 *            - uint32_t: modulo for the reduction 
 *            - uint32_t *out: output buffer
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - uint32_t *in1: first input buffer
 *            - size_t in1_msk_stride: stride between shares
 *            - size_t in1_data_stride: stride between masked bits
 *            - uint32_t *in2: second input buffer
 *            - size_t in2_msk_stride: stride between shares
 *            - size_t in2_data_stride: stride between masked bits
 **************************************************/
void secadd_modp(size_t nshares, size_t kbits, uint32_t q, uint32_t *out,
                 size_t out_msk_stride, size_t out_data_stride,
                 const uint32_t *in1, size_t in1_msk_stride,
                 size_t in1_data_stride, const uint32_t *in2,
                 size_t in2_msk_stride, size_t in2_data_stride) {

  uint32_t s[(kbits + 1) * nshares];
  uint32_t sp[(kbits + 1) * nshares];

  secadd(nshares, kbits, kbits + 1, s, 1, nshares, in1, in1_msk_stride,
         in1_data_stride, in2, in2_msk_stride, in2_data_stride);

  secadd_constant(nshares, kbits + 1, kbits + 1, sp, 1, nshares, s, 1, nshares,
                  (1 << (kbits + 1)) - q);

  secadd_constant_bmsk(nshares, kbits, kbits, out, out_msk_stride,
                       out_data_stride, sp, 1, nshares, q, &sp[kbits * nshares],
                       1);
}

/*************************************************
 * Name:        seca2b
 *
 * Description: Inplace arithmetic to boolean masking conversion:
 *            sum(in_i)%(2**kbits) = XOR(in_i) 
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words. Arithmetic
 *            masking is on 2**kbits.
 *            - uint32_t *in: input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 **************************************************/
void seca2b(size_t nshares, size_t kbits, uint32_t *in, size_t in_msk_stride,
            size_t in_data_stride) {


  size_t i, d;

  if (nshares == 1) {
    return;
  }

  size_t nshares_low = nshares / 2;
  size_t nshares_high = nshares - nshares_low;

  seca2b(nshares_low, kbits, in, in_msk_stride, in_data_stride);
  seca2b(nshares_high, kbits, &in[nshares_low * in_msk_stride], 
          in_msk_stride, in_data_stride);

  uint32_t expanded_low[kbits * nshares];
  uint32_t expanded_high[kbits * nshares];

  for (i = 0; i < kbits; i++) {
    for (d = 0; d < nshares_low; d++) {
      expanded_low[i * nshares + d] =
          in[i * in_data_stride + d * in_msk_stride];
      expanded_high[i * nshares + d] = 0;
    }
    for (d = nshares_low; d < nshares; d++) {
      expanded_high[i * nshares + d] =
          in[i * in_data_stride + d * in_msk_stride];
      expanded_low[i * nshares + d] = 0;
    }
  }

  secadd(nshares, kbits, kbits, in, in_msk_stride, in_data_stride, expanded_low,
         1, nshares, expanded_high, 1, nshares);
 
}

/*************************************************
 * Name:        seca2b_modp
 *
 * Description: Inplace arithmetic to boolean masking conversion:
 *            sum(in_i)%(p) = XOR(in_i) 
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words. kbits = ceil(log(p))
 *            - uint32_t p: modulus of the arithmetic masking
 *            - uint32_t *in: input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 **************************************************/
void seca2b_modp(size_t nshares, size_t kbits, uint32_t p, uint32_t *in,
                 size_t in_msk_stride, size_t in_data_stride) {

  size_t i, d;

  if (nshares == 1) {
    return;
  }

  size_t nshares_low = nshares / 2;
  size_t nshares_high = nshares - nshares_low;

  seca2b_modp(nshares_low, kbits, p, in, in_msk_stride, in_data_stride);
  seca2b_modp(nshares_high, kbits, p, &in[nshares_low * in_msk_stride], 
              in_msk_stride, in_data_stride);

  uint32_t expanded_low[(kbits + 1) * nshares];
  uint32_t expanded_high[(kbits + 1) * nshares];
  uint32_t u[(kbits + 1) * nshares];

  secadd_constant(nshares_low, kbits, kbits + 1, expanded_low, 1, nshares, in,
                  in_msk_stride, in_data_stride, (1 << (kbits + 1)) - p);

  for (i = 0; i < (kbits + 1); i++) {
    for (d = 0; d < nshares_low; d++) {
      // has already been written by secadd_constant_bmsk
      // expanded_low[i*nshares + d] = in[i*in_data_stride + d*in_msk_stride];
      expanded_high[i * nshares + d] = 0;
    }
    for (d = nshares_low; d < nshares; d++) {
      // kbits + 1 within in is unset
      expanded_high[i * nshares + d] =
          (i < (kbits)) ? in[i * in_data_stride + d * in_msk_stride] : 0;
      expanded_low[i * nshares + d] = 0;
    }
  }

  secadd(nshares, kbits + 1, kbits + 1, u, 1, nshares, expanded_high, 1,
         nshares, expanded_low, 1, nshares);

  secadd_constant_bmsk(nshares, kbits, kbits, in, in_msk_stride, in_data_stride,
                       u, 1, nshares, p, &u[kbits * nshares], 1);


}

/*************************************************
 * Name:        secb2a
 *
 * Description: Inplace boolean to arithmetic masking conversion:
 *            sum(in_i)%((1<<kbits)) = XOR(in_i) 
 *
 * /!\ This function is operating on two slices such that 64 conversions are
 * performed in parallel.
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words. kbits < 16
 *            - uint32_t p: modulus of the arithmetic masking
 *            - uint32_t *in: input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 **************************************************/
void secb2a(size_t nshares,
            size_t kbits,
            uint32_t *in, size_t in_msk_stride, size_t in_data_stride) {

  uint16_t z_dense[2 * BSSIZE * nshares];
  uint16_t zp_dense[2 * BSSIZE * nshares];
  uint32_t zp_str[2 * kbits * nshares];
  uint32_t b_str[2 * kbits * nshares];
  size_t d, i;
  uint32_t r;
  uint16_t r0, r1;
  uint32_t p = 1 << kbits;
  // generate uniform sharing for z
  // zp = p - z;
  for (d = 0; d < nshares - 1; d++) {
    for (i = 0; i < BSSIZE; i += 2) {
      r = get_random();
      r0 = r & ((1 << kbits) - 1);
      r1 = (r >> 16) & ((1 << kbits) - 1);
      z_dense[i * nshares + d] = r0;
      zp_dense[i * nshares + d] = p - r0;

      z_dense[(i + 1) * nshares + d] = r1;
      zp_dense[(i + 1) * nshares + d] = p - r1;

      r = get_random();
      r0 = r & ((1 << kbits) - 1);
      r1 = (r >> 16) & ((1 << kbits) - 1);
      z_dense[i * nshares + d + (BSSIZE * nshares)] = r0;
      zp_dense[i * nshares + d + (BSSIZE * nshares)] = p - r0;

      z_dense[(i + 1) * nshares + d + (BSSIZE * nshares)] = r1;
      zp_dense[(i + 1) * nshares + d + (BSSIZE * nshares)] = p - r1;
    }
#if ((BSSIZE & 0x1))
#error "Can only handle even number of slices."
#endif
  }

  // map zp to bitslice representation
  masked_dense2bitslice_opt(nshares - 1, kbits, zp_str, 1, nshares, (int16_t*)zp_dense, 1,
                            nshares);

  // last shares of zp set to zero
  for (i = 0; i < kbits; i++) {
    zp_str[i * nshares + (nshares - 1)] = 0;
    zp_str[i * nshares + (nshares - 1) + kbits * nshares] = 0;
  }

  // last shares of zp_str to zero
  seca2b(nshares, kbits, zp_str, 1, nshares);

  secadd(nshares, kbits, kbits, b_str, 1, nshares, in, in_msk_stride,
         in_data_stride, zp_str, 1, nshares);

  // map z to bistlice in output buffer
  masked_dense2bitslice_opt(nshares - 1, kbits, in, in_msk_stride,
                            in_data_stride, (int16_t*)z_dense, 1, nshares);

  // unmask b_str and set to the last share of the output
  for (i = 0; i < kbits; i++) {
    RefreshIOS_rec(nshares,nshares,&b_str[i * nshares], 1);

    in[i * in_data_stride + (nshares - 1) * in_msk_stride] = 0;
    for (d = 0; d < nshares; d++) {
      in[i * in_data_stride + (nshares - 1) * in_msk_stride] ^=
          b_str[i * nshares + d];
    }
  }
}

/*************************************************
 * Name:        RefreshIOS_rec
 *
 * Description: IOS refresh on boolean sharing 
 * 
 * Arguments: - size_t nshares: number of shares
 *            - size_t d: current recursion shars:
 *            - uint32_t *x: input buffer
 *            - size_t x_msk_stride: stride between shares
 **************************************************/
void RefreshIOS_rec(size_t nshares, size_t d,
    uint32_t *x, size_t x_msk_stride) {
  uint32_t r;
  if (d == 1) {
  } else if (d == 2) {
    r = get_random();
    x[0*x_msk_stride] ^= r;
    x[1*x_msk_stride] ^= r;
  } else {
    RefreshIOS_rec(nshares, d / 2, x, x_msk_stride);
    RefreshIOS_rec(nshares, d - d / 2, &x[(d / 2)*x_msk_stride], x_msk_stride);
    for (unsigned int i = 0; i < d / 2; i++) {
      r = get_random();
      x[i * x_msk_stride] ^= r;
      x[(i + d / 2) * x_msk_stride] ^= r;
    }
  }
}

/*************************************************
 * Name:        seccompress
 *
 * Description: Performs polynomial coefficients. See for details: 
 *              https://eprint.iacr.org/2021/1615.pdf, Algo 6.
 *
 * /!\ performs compression on 32 coefficients at once.
 *
 * Arguments: - size_t nshares: number of shares
 *            - uint32_t q: prime
 *            - uint32_t c: compression factor
 *            - uint32_t *out: output bitslice buffer. c-bits buffer.
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - int16_t *in: input buffer. Contains 32 coefficients
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 **************************************************/
void seccompress(size_t nshares, uint32_t q, uint32_t c,
                 uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
                 const int16_t *in, size_t in_msk_stride,
                 size_t in_data_stride) {

  size_t i, d;
  uint32_t ell = 0;
  uint32_t prod = q * nshares;
  while (prod > 0) {
    prod = prod >> 1;
    ell++;
  }

  uint32_t in_expanded[BSSIZE * nshares];
  uint32_t bs_expanded[(ell + c) * nshares];

  // map mod q to mod 2^ell.
  uint32_t tmp32;
  uint64_t tmp64;
  for (i = 0; i < BSSIZE; i++) {
    tmp64 = (in[i * in_data_stride] + q) % q;
    tmp32 = ((((tmp64 << (c + ell + 1)) + q) >> 1) / (q));
    tmp32 += (1 << (ell - 1));
    in_expanded[i * nshares] = tmp32 & ((1 << (ell + c)) - 1);
    for (d = 1; d < nshares; d++) {
      tmp64 = (q + in[i * in_data_stride + d * in_msk_stride]) % q;
      tmp32 = ((((tmp64 << (c + ell + 1)) + q) >> 1) / (q));
      in_expanded[i * nshares + d] = tmp32 & ((1 << (ell + c)) - 1);
    }
  }

  // map to bitslice
  masked_dense2bitslice_opt_u32(nshares, ell + c, bs_expanded, 1, nshares,
                            in_expanded, 1, nshares);
  
  // convert A2B
  seca2b(nshares, ell + c, bs_expanded, 1, nshares);

  // map to the output
  for (i = 0; i < c; i++) {
    for (d = 0; d < nshares; d++) {
      out[i * out_data_stride + d * out_msk_stride] =
          bs_expanded[(ell + i) * nshares + d];
    }
  }
}

/*************************************************
 * Name:        secfulladd
 *
 * Description: Performs a 3-bit full adder. See for details: 
 *              https://eprint.iacr.org/2022/158.pdf, Algo 5.
 *
 * Arguments: - size_t nshares: number of shares
 *            - uint32_t *co: carry out buffer
 *            - size_t co_msk_stride: carry out shares stride
 *            - uint32_t *w: output bit
 *            - size_t w_msk_stride: w shares stride
 *            - uint32_t *ci: carry in buffer
 *            - size_t ci_msk_stride: carry in shares stride
 *            - uint32_t *x: first input bit
 *            - size_t x_msk_stride: x shares stride
 *            - uint32_t *y: second input bit
 *            - size_t y_msk_stride: y shares stride
 **************************************************/
void secfulladd(size_t nshares, uint32_t *co, size_t co_msk_stride, uint32_t *w,
                size_t w_msk_stide, const uint32_t *ci, size_t ci_msk_stride,
                const uint32_t *x, size_t x_msk_stride, const uint32_t *y,
                size_t y_msk_stride) {
  uint32_t a[nshares];
  uint32_t b[nshares];

  masked_xor(nshares, a, 1, x, x_msk_stride, y, y_msk_stride);

  masked_xor(nshares, b, 1, x, x_msk_stride, ci, ci_msk_stride);

  // compute carry out
  masked_and(nshares, a, 1, b, 1, a, 1);

  masked_xor(nshares, a, 1, a, 1, x, x_msk_stride);

  // compute w
  masked_xor(nshares, w, w_msk_stide, b, 1, y, y_msk_stride);

  // write carry out
  for (size_t d = 0; d < nshares; d++) {
    co[d * co_msk_stride] = a[d];
  }
}

/*************************************************
 * Name:        sechalfadd
 *
 * Description: Performs a 2-bit half adder.
 *
 * Arguments: - size_t nshares: number of shares
 *            - uint32_t *c: carry out buffer
 *            - size_t c_msk_stride: carry out shares stride
 *            - uint32_t *w: output bit
 *            - size_t w_msk_stride: w shares stride
 *            - uint32_t *x: first input bit
 *            - size_t x_msk_stride: x shares stride
 *            - uint32_t *y: second input bit
 *            - size_t y_msk_stride: y shares stride
 **************************************************/
void sechalfadd(size_t nshares, uint32_t *c, size_t c_msk_stride, uint32_t *w,
                size_t w_msk_stide,
                const uint32_t *x, size_t x_msk_stride, const uint32_t *y,
                size_t y_msk_stride) {
  uint32_t a[nshares];

  masked_xor(nshares, a, 1, x, x_msk_stride, y, y_msk_stride);

  masked_and(nshares, c, c_msk_stride, x, x_msk_stride, y, y_msk_stride);
  
  // write result out
  for (size_t d = 0; d < nshares; d++) {
    w[d * w_msk_stide] = a[d];
  }
}

/*************************************************
 * Name:        sechalfsub
 *
 * Description: Performs a 2-bit half subtracter.
 *
 * Arguments: - size_t nshares: number of shares
 *            - uint32_t *c: carry out buffer
 *            - size_t c_msk_stride: carry out shares stride
 *            - uint32_t *w: output bit
 *            - size_t w_msk_stride: w shares stride
 *            - uint32_t *x: first input bit
 *            - size_t x_msk_stride: x shares stride
 *            - uint32_t *y: second input bit
 *            - size_t y_msk_stride: y shares stride
 **************************************************/
void sechalfsub(size_t nshares, uint32_t *c, size_t c_msk_stride, uint32_t *w,
                size_t w_msk_stide,
                const uint32_t *x, size_t x_msk_stride, const uint32_t *y,
                size_t y_msk_stride) {
  uint32_t a[nshares], notx[nshares];
  
  notx[0] = ~x[0];
  for (size_t i = 1; i < nshares; i++)
  {
    notx[i] = x[i * x_msk_stride];
  }

  masked_xor(nshares, a, 1, x, x_msk_stride, y, y_msk_stride);

  masked_and(nshares, c, c_msk_stride, notx, 1, y, y_msk_stride);
  
  // write result out
  for (size_t d = 0; d < nshares; d++) {
    w[d * w_msk_stide] = a[d];
  }
}


/*************************************************
 * Name:        secb2a_32to48
 *
 * Description: Inplace boolean to arithmetic masking conversion:
 *            sum(in_i)%((1<<64)) = XOR(in_i) 
 *
 * Arguments: - size_t nshares: number of shares
 *            - uint32_t *in: input buffer (48 values of width 32 bit, where the lower 32 values are the input, 
 *                 upper 16 values can be initialized arbitrarily and are overwritten)
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 **************************************************/
void secb2a_32to48(size_t nshares, uint32_t *in, size_t in_msk_stride, size_t in_data_stride)
{
  uint32_t z[nshares-1][64]; // we still work on 64-bit arrays here to enable usage of transpose32
  uint32_t zp[nshares][64];
  uint32_t b_str[48 * nshares];
  size_t d, i;
  uint64_t r,rp;
  // generate uniform sharing for z
  // zp = (1<<64)-z;
  for (d = 0; d < nshares - 1; d++)
  {
    for (i = 0; i < 32; i += 2)
    {
      r = (((uint64_t)get_random())<<32) | get_random();
      rp = ((uint64_t)1<<48)-(r&(((uint64_t)1<<48)-1));
      z[d][i] = r&(((uint64_t)1<<48)-1);
      z[d][i + 32] = (r&(((uint64_t)1<<48)-1))>>32;
      zp[d][i] = rp;
      zp[d][i + 32] = rp>>32;

      r &= (uint64_t)0xffff000000000000; // cancel lower part
      r >>= 16; // upper 16 bit of r are still fresh randomness
      r |= get_random(); // add fresh randomness for lower part
      rp = ((uint64_t)1<<48)-r;
      z[d][i + 1] = r;
      z[d][i + 1 + 32] = r>>32;
      zp[d][i + 1] = rp;
      zp[d][i + 1 + 32] = rp>>32;
    }
    transpose32(z[d]);
    transpose32(&z[d][32]);
    transpose32(zp[d]);
    transpose32(&zp[d][32]);
  }

  // last shares of zp set to zero
  for (i = 0; i < 48; i++) 
  {
    zp[nshares-1][i] = 0;
  }

  // last shares of zp_str to zero
  seca2b(nshares, 48, &zp[0][0], 64, 1);
  secadd(nshares, 32, 33, b_str, 1, nshares, in, in_msk_stride, in_data_stride, &zp[0][0], 64, 1);
  // carry now is in b_str[32 * nshares], we need to add it to zp
  for (i = 32; i < 47; i++)
  { 
    sechalfadd(nshares, //size_t nshares, 
               &b_str[(i+1) * nshares], //uint32_t *c,
               1, //size_t c_msk_stride, 
               &b_str[i * nshares], //uint32_t *w,
               1, //size_t w_msk_stide,
               &b_str[i * nshares], //uint32_t *x,
               1, //size_t x_msk_stide,
               &zp[0][i], //const uint32_t *y, 
               64 //size_t y_msk_stride
              );
  }
  masked_xor(nshares, &b_str[47 * nshares], 1, &b_str[47 * nshares], 1, &zp[0][47], 64);
   
   
  // map z to bistlice in output buffer
  for (d = 0; d < nshares - 1; d++)
  { 
    for (i = 0; i < 48; i++) 
    {
      in[i * in_data_stride + d * in_msk_stride] = z[d][i];
    }
  }
  
  // unmask b_str and set to the last share of the output
  for (i = 0; i < 48; i++)
  {
    RefreshIOS_rec(nshares, nshares, &b_str[i * nshares], 1);
      
    in[i * in_data_stride + (nshares - 1) * in_msk_stride] = 0;
    for (d = 0; d < nshares; d++)
    {
      in[i * in_data_stride + (nshares - 1) * in_msk_stride] ^= b_str[i * nshares + d];
    }
  }
}


// old = cond ? new : old
// expects cond to be one bit
void masked_sel(uint32_t *cond, size_t cond_stride, uint32_t *old, size_t old_stride, 
                uint32_t *new, size_t new_stride) {
  uint32_t cond_mask[NSHARES];
  uint32_t xor[NSHARES];

  for(int s=0; s<NSHARES; s++) {
    cond_mask[s] = ~(cond[s*cond_stride] -1); // expand from 1 bit to registerwidth
  }
  // xor = old ^ new -> if cond : old = old ^ old ^ new = new, else : old = old ^ 0 = old
  masked_xor(NSHARES, xor, 1, old, old_stride, new, new_stride);
  masked_and(NSHARES, xor, 1, xor, 1, cond_mask, 1);
  masked_xor(NSHARES, old, old_stride, old, old_stride, xor, 1);
}

// old = cond ? new : old
// bits <= 32
// expects datastride of 1
void masked_bs_sel(uint32_t *cond, size_t cond_stride, uint32_t *old, size_t old_stride, 
                uint32_t *new, size_t new_stride, size_t bits) {
  uint32_t xor[NSHARES];
  
  for(int i=0; i<bits; i++) {
    masked_xor(NSHARES, xor, 1, &old[i], old_stride, &new[i], new_stride);
    masked_and(NSHARES, xor, 1, xor, 1, cond, cond_stride);
    masked_xor(NSHARES, &old[i], old_stride, &old[i], old_stride, xor, 1);    
  }
}

// bits <= 32
// expects datastride of 1
// returns 1 per element if equal
void masked_bs_eq(uint32_t *res, size_t res_stride, uint32_t *a, size_t a_stride,
                  uint32_t *b, size_t b_stride, size_t bits) {
  uint32_t xor[NSHARES][bits];

  //xor = a^b: (if a==b: xor==0, else xor != 0)
  for(int s=0; s<NSHARES; s++) {
    for(int i=0; i<bits; i++) {
      xor[s][i] = a[s*a_stride + i] ^ b[s*b_stride + i];
    }
  }

  // OR of all bits in xor is 0 if xor==0, and is 1 otherwise
  for(int s=0; s<NSHARES; s++) {
    res[s*res_stride] = 0; // init res with 0
  }

  for(int i=0; i<bits; i++) {
    masked_or(NSHARES, res, res_stride, res, res_stride, &xor[0][i], bits);
  }
  res[0] = ~res[0]; //flip the result: if equal : xor = 0 -> res = 0 -> res = 1
}

// returns 1 if b > a
void cmpg(uint32_t *res, const uint32_t *a, uint32_t stride_a, const uint32_t *b, uint32_t stride_b, size_t bits)
{
  size_t n,i;
  uint32_t B[NSHARES], xor[NSHARES], tmp[NSHARES];
  
  masked_xor(NSHARES, xor, 1, a, stride_a, b, stride_b);
  
  for (n = 0; n < NSHARES; n++)
  {
    B[n] = b[n*stride_b];
    res[n] = 0;
    tmp[n] = 0;
  }
  
  for (i = 0; i < bits; i++) 
  {
    masked_xor(NSHARES, tmp, 1, res, 1, B, 1); // tmp_i = res ^ b_i
    masked_and_1bit(NSHARES, tmp, 1, tmp, 1, xor, 1); // tmp_i = tmp_i & xor_i
    masked_xor(NSHARES, res, 1, res, 1, tmp, 1); // res = res ^ tmp_i = res ^ ((res ^ b_i) & xor_i)
    
    for (n = 0; n < NSHARES; n++)
    {
      xor[n] >>= 1;
      B[n] >>= 1;
    }
  }

  for (n = 0; n < NSHARES; n++)
  {
    res[n] = res[n] & 1; // return only lowest bit
  }
}

// returns 1 if b > a
void bscmpg(uint32_t *res, const uint32_t *a, const uint32_t stride_a, const uint32_t *b, const uint32_t stride_b)
{
  size_t n,i;
  uint32_t xor[NSHARES], tmp[NSHARES];
  
  for (n = 0; n < NSHARES; n++)
  {
    res[n] = 0;
  }
  
  for (i = 0; i < LOG_N; i++) 
  {
    masked_xor(NSHARES, xor, 1, a + i, stride_a, b + i, stride_b);
    masked_xor(NSHARES, tmp, 1, res, 1, b + i, stride_b);
    masked_and(NSHARES, tmp, 1, tmp, 1, xor, 1);
    masked_xor(NSHARES, res, 1, res, 1, tmp, 1);
  }
}
