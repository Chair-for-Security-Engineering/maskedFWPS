#include <stdint.h>
#include <stddef.h>
#include "gadgets.h"
#include "params.h"
#include "masked_fwsampling.h"

#if N == 12323 // BIKE-1
#define L 9
#elif N == 24659 // BIKE-3
#define L 8
#elif N == 40973 // BIKE-5
#define L 11
#elif N == 17669 // HQC-128
#define L 8
#elif N == 35851 // HQC-192
#define L 10
#elif N == 57637 // HQC-256
#define L 12
#elif N == 3488 // McE 348864
#define L 8
#elif N == 4608 // McE 460896
#define L 9
#elif N == 6688 // McE 6688128
#define L 8
#elif N == 6960 // McE 6960119
#define L 9
#elif N == 8192 // McE 8192128
#define L 6
#elif N == 509
#define L 1
#elif N == 677
#define L 3
#elif N == 821
#define L 3
#elif N == 653
#define L 4
#elif N == 761
#define L 3
#elif N == 857
#define L 3
#elif N == 953
#define L 5
#elif N == 1013
#define L 4
#elif N == 1277
#define L 3
#else
#error
#endif

#ifdef NTRULPR
#if N == 653
#define L 3
#elif N == 761
#define L 6
#elif N == 857
#define L 3
#elif N == 953
#define L 3
#elif N == 1013
#define L 3
#elif N == 1277
#define L 5
#else
#error
#endif
#endif

static void sample_fixed_32(uint32_t *out, size_t stride, uint32_t randomness[NSHARES][L])
{
#ifdef NTRULPR
#if N == 653 || N == 857 || N == 953 || N == 1013 // t = 3
  masked_and(NSHARES, out, stride, &randomness[0][0], L, &randomness[0][1], L);
  out[0] = ~out[0];
  for (size_t i = 2; i < L; i++)
  {
    masked_and(NSHARES, out, stride, &randomness[0][i], L, out, stride);
  }
#elif N == 761 // t = 21
  masked_and(NSHARES, out, stride, &randomness[0][0], L, &randomness[0][1], L);
  out[0] = ~out[0];
  masked_and(NSHARES, out, stride, &randomness[0][2], L, out, stride);
  out[0] = ~out[0];
  masked_and(NSHARES, out, stride, &randomness[0][3], L, out, stride);
  out[0] = ~out[0];
  masked_and(NSHARES, out, stride, &randomness[0][4], L, out, stride);
  out[0] = ~out[0];
  masked_and(NSHARES, out, stride, &randomness[0][5], L, out, stride);
#elif N == 1277 // t = 11
  masked_and(NSHARES, out, stride, &randomness[0][0], L, &randomness[0][1], L);
  out[0] = ~out[0];
  masked_and(NSHARES, out, stride, &randomness[0][2], L, out, stride);
  out[0] = ~out[0];
  masked_and(NSHARES, out, stride, &randomness[0][3], L, out, stride);
  out[0] = ~out[0];
  for (size_t i = 4; i < L; i++)
  {
    masked_and(NSHARES, out, stride, &randomness[0][i], L, out, stride);
  }
#else
#error
#endif
#else
#if N == 12323 || N == 35851 || N == 677 || N == 761 || N == 857 || N == 1277 // t = 3
  masked_and(NSHARES, out, stride, &randomness[0][0], L, &randomness[0][1], L);
  out[0] = ~out[0];
  for (size_t i = 2; i < L; i++)
  {
    masked_and(NSHARES, out, stride, &randomness[0][i], L, out, stride);
  }
#elif N == 24659 || N == 17669 || N == 8192 // t = 1
  masked_and(NSHARES, out, stride, &randomness[0][0], L, &randomness[0][1], L);
  for (size_t i = 2; i < L; i++)
  {
    masked_and(NSHARES, out, stride, &randomness[0][i], L, out, stride);
  }
#elif N == 40973 || N == 1013 || N == 653 // t = 7
  masked_and(NSHARES, out, stride, &randomness[0][0], L, &randomness[0][1], L);
  masked_and(NSHARES, out, stride, &randomness[0][2], L, out, stride);
  out[0] = ~out[0];
  for (size_t i = 3; i < L; i++)
  {
    masked_and(NSHARES, out, stride, &randomness[0][i], L, out, stride);
  }
#elif N == 57637 || N == 6960 // t = 9
  masked_and(NSHARES, out, stride, &randomness[0][0], L, &randomness[0][1], L);
  masked_and(NSHARES, out, stride, &randomness[0][2], L, out, stride);
  out[0] = ~out[0];
  masked_and(NSHARES, out, stride, &randomness[0][3], L, out, stride);
  out[0] = ~out[0];
  for (size_t i = 4; i < L; i++)
  {
    masked_and(NSHARES, out, stride, &randomness[0][i], L, out, stride);
  }
#elif N == 3488 || N == 6688 || N == 821 // t = 5
  masked_and(NSHARES, out, stride, &randomness[0][0], L, &randomness[0][1], L);
  out[0] = ~out[0];
  masked_and(NSHARES, out, stride, &randomness[0][2], L, out, stride);
  out[0] = ~out[0];
  for (size_t i = 3; i < L; i++)
  {
    masked_and(NSHARES, out, stride, &randomness[0][i], L, out, stride);
  }
#elif N == 4608 // t = 11
  masked_and(NSHARES, out, stride, &randomness[0][0], L, &randomness[0][1], L);
  out[0] = ~out[0];
  masked_and(NSHARES, out, stride, &randomness[0][2], L, out, stride);
  out[0] = ~out[0];
  masked_and(NSHARES, out, stride, &randomness[0][3], L, out, stride);
  out[0] = ~out[0];
  for (size_t i = 4; i < L; i++)
  {
    masked_and(NSHARES, out, stride, &randomness[0][i], L, out, stride);
  }
#elif N == 509 // t = 1
  for (size_t d = 0; d < NSHARES; d++)
  {
    out[d*stride] = randomness[d][0];
  }
#elif N == 953
  masked_and(NSHARES, out, stride, &randomness[0][0], L, &randomness[0][1], L);
  out[0] = ~out[0];
  masked_and(NSHARES, out, stride, &randomness[0][2], L, out, stride);
  masked_and(NSHARES, out, stride, &randomness[0][3], L, out, stride);
  out[0] = ~out[0];
  masked_and(NSHARES, out, stride, &randomness[0][4], L, out, stride);
#else
#error
#endif
#endif
}

uint32_t inputrandbuffer[4099];
uint32_t inputrand(void)
{
  static size_t offset = 0;
  uint32_t rv = inputrandbuffer[offset];
  offset += 1;
  offset %= 4099;
  return rv;
}

// fixed weight sampling for bike
uint8_t sample_fixed(uint32_t out[NSHARES][(N+31)/32]) // return weight to prevent optimization
{
  uint32_t r[NSHARES][L];
  uint32_t out_arith_bs[NSHARES][8], out_arith_bs_x4[32], negA0[32];
  size_t bsoffset = 0;
  uint8_t weight = 0;
  //while (weight != 71)
  {
    uint8_t out_weight[NSHARES] = {0};
    uint16_t coefcnt[NSHARES] = {0};
    // sample polynomial with each coefficient having probability 3/512 of being one
    // which approximates the correct value 71/12323
    for (size_t i = 0; i < (N+31)/32; i++)
    {
      for (size_t j = 0; j < NSHARES*L; j++)
      {
        r[j/L][j%L] = inputrand();
      }
      sample_fixed_32(&out[0][i], (N+31)/32, r);
    }
    
    // convert to arithmetic shares mod 256 and accumulate weight
    for (size_t i = 0; i < (N+31)/32; i++)
    {
      // B2A conversion
      uint32_t B1[2], carry[NSHARES];

      B1[1] = get_random();
      negA0[bsoffset+0] = get_random();
      B1[0] = negA0[bsoffset+0]^B1[1];
      sechalfadd(NSHARES,
                  carry, 1,
                  &out_arith_bs[0][0], 8,
                  &out[0][i], (N+31)/32,
                  B1, 1);
      for (size_t j = 1; j < 7; j++)
      {
        B1[1] = get_random();
        negA0[bsoffset+j] = get_random();
        B1[0] = negA0[bsoffset+j]^B1[1];
        sechalfadd(NSHARES,
                   carry, 1,
                   &out_arith_bs[0][j], 8,
                   carry, 1,
                   B1, 1);
      }
      B1[1] = get_random();
      negA0[bsoffset+7] = get_random();
      B1[0] = negA0[bsoffset+7]^B1[1];
      masked_xor(NSHARES,
                 &out_arith_bs[0][7], 8,
                 carry, 1,
                 B1, 1);

      // now: each slice in out_arith_bs[0]^out_arith_bs[1] is one arithmetic share; 256 minus each slice in negA0 is the other share
      for (size_t j = 0; j < 8; j++)
      {
        out_arith_bs_x4[bsoffset+j] = out_arith_bs[0][j];
        for (size_t n = 1; n < NSHARES; n++)
        {
          out_arith_bs_x4[bsoffset+j] ^= out_arith_bs[n][j];
        }
      }

      bsoffset += 8;
      if (bsoffset == 32 || i == ((N+31)/32-1))
      {
        // un-bitslice
        bsoffset = 0;

        transpose32(out_arith_bs_x4);
        for (size_t j = 0; j < 4; j++)
        {
          for (size_t k = 0; k < 32; k++)
          {
            if (coefcnt[0] == N)
              break;
            out_weight[0] += (out_arith_bs_x4[k] >> (j*8))&0xff;
            coefcnt[0] += 1;
          }
        }

        transpose32(negA0);
        for (size_t j = 0; j < 4; j++)
        {
          for (size_t k = 0; k < 32; k++)
          {
            if (coefcnt[1] == N)
              break;
            out_weight[1] -= (negA0[k] >> (j*8))&0xff;
            coefcnt[1] += 1;
          }
        }
      }
    }

    // unmask accumulated weight and store to `weight` variable, then this is repeated until the correct weight is found
    weight = (out_weight[0] + out_weight[1]) % 256;
    bsoffset = 0;
  }
}
