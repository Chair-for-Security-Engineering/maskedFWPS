#include <stdint.h>
#include <stddef.h>
#include "gadgets.h"
#include "params.h"
#include "masked_fwsampling.h"

//expects randomness in p for first W entries
void bs_fisher_yates(uint32_t p[NSHARES][W_PADDED]) {
  
  uint32_t t[NSHARES][64];
  union U64 x;

  for(int s=0; s<NSHARES; s++) {
    for(int i=0; i<64; i++) {
      t[s][i] = 0; // init t
    }
  }

  //component 1: sampling int out of range: p[i] = i + rand(n-1-i);
  for(int i=0; i<W_PADDED; i+=32) {
    //assume we receive our random values boolean bitsliced masking domain 
    for(int s=0; s<NSHARES; s++) {
      for(int j=0; j<32; j++) {
        t[s][j] = p[s][i+j]; // TODO: make b2a non-inplace to get rid of copying
      }
    }
    secb2a_32to48(NSHARES, &t[0][0], 64, 1); // transform to mod 2^48 arithmetic domain;
    transpose32_array(&t[0][0], 64); // unbitslice

    for(int j=0; j<32; j++) {
      for(int s=0; s<NSHARES; s++) {
        x.v32[0] = t[s][j+ 0];  // lower 32bit
        x.v32[1] = t[s][j+32];  // upper 32bit
        x.v64 = x.v64 * ((uint64_t) N-1-i-j);
        t[s][j+ 0] = x.v32[0];
        t[s][j+32] = x.v32[1] & 0x0000FFFF; //mod 2^48
      }
      t[0][j+32] = (t[0][j+32] + i+j) & 0x0000FFFF; // add i to the upper half, compute mod 2^48
    }
    //a2b:
    transpose32_array(&t[0][0], 64); // bitslice
    seca2b(NSHARES, 48, &t[0][0], 64, 1);
    for(int s=0; s<NSHARES; s++) {
      for(int j=0; j<32; j++) {
        p[s][i+j] = t[s][j+32]; // copy back upper 32bits
      }
    }
  }

  //component 2: if collision then exchange bitsliced
  int start_j32 = W_PADDED-32; // init with rightmost index in 32er steps
  int i32 = W_PADDED-32;

  uint32_t bs_i[NSHARES][32]; // 32 times i
  uint32_t bs_pi[NSHARES][32]; // 32 times p[i]
  uint32_t res[NSHARES]; // result of comparison

  for(int s=1; s<NSHARES; s++) { // i is public, no randomness in other shares
    for(int j=0; j<32; j++) {
        bs_i[s][j] = 0;
    }
    // bitslicing zeros would not change anything
  }

  for(int i=W-2; i>=0; i--) {
    start_j32 = start_j32 > i+1 ? start_j32-32 : start_j32; // index of 32er block containing j
    i32 = i32 > i ? i32-32 : i32; // index of 32er block containing i

    // prepare 32times index i:
    for(int j=0; j<32; j++) {
      bs_i[0][j] = i;
    }
    transpose32(&bs_i[0][0]);

    // prepare 32times p[i]:
    for(int s=0; s<NSHARES; s++) {
      transpose32(&p[s][i32]); // unbitslice to avoid bitfiddeling
      for(int j=0; j<32; j++) {
        bs_pi[s][j] = p[s][i]; // copy p[i] 32 times for comparison
      }
      transpose32(&p[s][i32]); // bitslice again
      transpose32(&bs_pi[s][0]);
    }

    for(int j32 = start_j32; j32+31 < W_PADDED; j32 += 32 ) {
      int j=i+1;
      int diff = j-j32;
      uint32_t mask = 0xFFFFFFFF;
      if(diff > 0)
        mask = mask << diff; //rightmost diff bits are zero -> dont touch these indizes

      masked_bs_eq(res, 1, &p[0][j32], W_PADDED, &bs_pi[0][0], 32, LOG_N);
      for(int s=0; s<NSHARES; s++) {
        res[s] &= mask; // indizes we should not touch in this iteration are set to 0
      }
      masked_bs_sel(res, 1, &p[0][j32], W_PADDED, &bs_i[0][0], 32, LOG_N);
    }
  }

  for(int i=0; i<W_PADDED; i+=32) {
    transpose32_shares(&p[0][i], W_PADDED);
  }
}
