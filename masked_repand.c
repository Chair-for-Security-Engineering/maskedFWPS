#include <stdint.h>
#include <stddef.h>
#include "gadgets.h"
#include "params.h"
#include "masked_fwsampling.h"

volatile uint32_t numwcheck = 0;
static uint16_t masked_weight_check(const uint32_t out[NSHARES][(N+31)/32])
{
  numwcheck += 1;
  uint16_t out_weight[NSHARES] = {0};
  uint16_t coefcnt[NSHARES] = {0};
  uint32_t out_arith_bs[NSHARES][16], out_arith_bs_x3[32], negA0[32];
  size_t bsoffset = 0;
  // convert to arithmetic shares mod 65536 and accumulate weight
  for (size_t i = 0; i < (N+31)/32; i++)
  {
    // B2A conversion
    uint32_t B1[2], carry[NSHARES];

    B1[1] = get_random();
    negA0[bsoffset+0] = get_random();
    B1[0] = negA0[bsoffset+0]^B1[1];
    sechalfadd(NSHARES,
                carry, 1,
                &out_arith_bs[0][0], 16,
                &out[0][i], (N+31)/32,
                B1, 1);
    for (size_t j = 1; j < 11; j++)
    {
      B1[1] = get_random();
      negA0[bsoffset+j] = get_random();
      B1[0] = negA0[bsoffset+j]^B1[1];
      sechalfadd(NSHARES,
                 carry, 1,
                 &out_arith_bs[0][j], 16,
                 carry, 1,
                 B1, 1);
    }
    B1[1] = get_random();
    negA0[bsoffset+11] = get_random();
    B1[0] = negA0[bsoffset+11]^B1[1];
    masked_xor(NSHARES,
               &out_arith_bs[0][11], 16,
               carry, 1,
               B1, 1);

    // now: each slice in out_arith_bs[0]^out_arith_bs[1] is one arithmetic share; 65536 minus each slice in negA0 is the other share
    for (size_t j = 0; j < 16; j++)
    {
      out_arith_bs_x3[bsoffset+j] = out_arith_bs[0][j];
      for (size_t n = 1; n < NSHARES; n++)
      {
        out_arith_bs_x3[bsoffset+j] ^= out_arith_bs[n][j];
      }
    }

    bsoffset += 16;
    if (bsoffset == 32 || i == ((N+31)/32-1))
    {
      // un-bitslice
      bsoffset = 0;

      transpose32(out_arith_bs_x3);
      for (size_t j = 0; j < 2; j++)
      {
        for (size_t k = 0; k < 32; k++)
        {
          if (coefcnt[0] == N)
            break;
          out_weight[0] += (out_arith_bs_x3[k] >> (j*16))&0xffff;
          coefcnt[0] += 1;
        }
      }

      transpose32(negA0);
      for (size_t j = 0; j < 2; j++)
      {
        for (size_t k = 0; k < 32; k++)
        {
          if (coefcnt[1] == N)
            break;
          out_weight[1] -= (negA0[k] >> (j*16))&0xffff;
          coefcnt[1] += 1;
        }
      }
    }
  }

  // unmask accumulated weight and store to `weight` variable, then this is repeated until the correct weight is found
  return (out_weight[0] + out_weight[1]) % (1<<12);
}

void repand_sample_fixed(uint32_t A[NSHARES][(N+31)/32], const uint16_t target_weight) // we assume A to be initialized to zero!
{
  size_t i,d;
  uint16_t w = target_weight;
  uint32_t Abar[NSHARES][(N+31)/32], Aj[NSHARES];

  
  do {
    for (i = 0; i < (N+31)/32; i++)
    {
      for (d = 0; d < NSHARES; d++)
      {
        Aj[d] = randomint();
      }
      A[0][i] = ~A[0][i];
      masked_and(NSHARES, &Abar[0][i], (N+31)/32, Aj, 1, &A[0][i], (N+31)/32);
      A[0][i] = ~A[0][i];
    }
    do {
      for (i = 0; i < (N+31)/32; i++)
      {
        for (d = 0; d < NSHARES; d++)
        {
          Aj[d] = randomint();
        }
        masked_and(NSHARES, &Abar[0][i], (N+31)/32, &Abar[0][i], (N+31)/32, Aj, 1);
      }
    } while (masked_weight_check(Abar) > w);
    for (i = 0; i < (N+31)/32; i++)
    {
      A[0][i] = ~A[0][i];
      Abar[0][i] = ~Abar[0][i];
      masked_and(NSHARES, &A[0][i], (N+31)/32, &A[0][i], (N+31)/32, &Abar[0][i], (N+31)/32);
      A[0][i] = ~A[0][i];
    }
    w = target_weight - masked_weight_check(A);
  } while(w != 0);
}
