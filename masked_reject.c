#include <stdint.h>
#include <stddef.h>
#include "gadgets.h"
#include "params.h"
#include "masked_fwsampling.h"

/////////////////////////////////////////////////////////////////////////////////////////////
// Bounded Rejection ////////////////////////////////////////////////////////////////////////

//expects non bitsliced input
//expects rand filled with randomness
int bs_ct_reject(uint32_t p[NSHARES][W_PADDED], uint32_t rand[NSHARES][BOUND]) {
  uint32_t cntr[NSHARES][LOG_BOUND];
  uint32_t dup[NSHARES];
  uint32_t tmp[NSHARES];
  uint32_t r[NSHARES][32];
  uint32_t f[NSHARES][32];
  uint32_t c32[NSHARES][32];
  uint32_t n_masked[NSHARES];

  // init arrays with value that does not create false collision:
  for(int s=0; s<NSHARES; s++) {
    for(int i=0; i<LOG_BOUND; i++) {
      cntr[s][i] = N;
    }
  }

  for(int s=0; s<NSHARES; s++) {
    for(int i=0; i<W_PADDED; i++) {
      p[s][i] = 0;
    }
  }

  n_masked[0] = N;
  for(int s=1; s<NSHARES; s++) {
    n_masked[s] = 0; // public value, zero masks;
  }

  for(int i=0; i<BOUND; i++) {
    // init dup and r:
    for(int s=0; s<NSHARES; s++) {
      dup[s] = 0;
      for(int j=0; j<32; j++) {
        r[s][j] = rand[s][i];
      } 
    }
    transpose32_shares(&r[0][0], 32);
    // compare and cmov bitsliced:
    for(int c=0; c<W_PADDED; c+=32) {
      // init c:
      for(int j=0; j<32; j++) {
        c32[0][j] = j+c;
      }
      for(int s=1; s<NSHARES; s++) {
        for(int j=0; j<32; j++) {
          c32[s][j] = 0; // public value, zero masks;
        }       
      }
      transpose32(&c32[0][0]);

      masked_bs_eq(tmp, 1, &r[0][0], 32, &p[0][c], W_PADDED, LOG_N+1); // tmp = p[c] == r 
      masked_or(NSHARES, dup, 1, dup, 1, tmp, 1); // dup |= tmp
      masked_bs_eq(tmp, 1, &c32[0][0], 32, &cntr[0][0], LOG_BOUND, LOG_W); // tmp = c == cntr
      masked_bs_sel(tmp, 1, &p[0][c], W_PADDED, &r[0][0], 32, LOG_N); // cmov(p[c], r, tmp)
    }

    //dup == 0 ?
    for(int j=1; j<32; j++) {
      for(int s=0; s<NSHARES; s++) {
        tmp[s] = (dup[s]>>j) &1;
      }
      masked_or_1bit(NSHARES, dup, 1, dup, 1, tmp, 1);
    }
    dup[0] = ~dup[0] &1; //dup = 1 if dup == 0

    cmpg(tmp, &rand[0][i], BOUND, n_masked, 1, LOG_N); // tmp = r < N
    masked_and(NSHARES, tmp, 1, tmp, 1, dup, 1); // tmp = dup == 0 & r<N
    
    for(int s=0; s<NSHARES; s++) {
      tmp[s] = ~((tmp[s]&1) + 0xffffffff); // extend 1 bit to 32
    }

    // cntr += t
    for(int j=0; j<LOG_BOUND; j++) {
      sechalfadd(NSHARES, tmp, 1, &cntr[0][j], LOG_BOUND, &cntr[0][j], LOG_BOUND, tmp, 1);
    }
  }
  transpose32_array(&p[0][0], W_PADDED); // unbitslice
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Simple Rejection /////////////////////////////////////////////////////////////////////////

//expects rand filled with randomness
int bs_reject(uint32_t p[NSHARES][W_PADDED], uint32_t rand[NSHARES][B_PADDED]) {
  uint32_t i = 0;
  uint32_t c = 0;
  uint32_t r = 0;
  uint32_t index;
  uint32_t r32[NSHARES][LOG_N];
  uint32_t n32[NSHARES][32];
  uint32_t valid[NSHARES];
  uint32_t res[NSHARES];

  //init N:
  for(int j=0; j<32; j++) {
    n32[0][j] = N;
  }
  for(int s=1; s<NSHARES; s++) {
    for(int j=0; j<32; j++) {
      n32[s][j] = 0;  // public value, zero masking
    }
  }
  transpose32(&n32[0][0]);

  //init p:
  for(int s=0; s<NSHARES; s++) {
    for(int j=0; j<W_PADDED; j++) {
      p[s][j] = 0;  // public value, zero masking
    }
  }
  transpose32_array(&p[0][0], W_PADDED);

  //init res:
  for(int s=0; s<NSHARES; s++) {
    res[s] = 0;
  }

  while(i < W) {
    // compare next 32 random values
    if(r>= B_PADDED)
      return 1; // very unlikely event of not enough randomness
    transpose32_shares(&rand[0][r], B_PADDED);
    bscmpg(valid, &rand[0][r], B_PADDED, &n32[0][0], 32);
    for(int s=1; s<NSHARES; s++) {
      valid[0] ^= valid[s]; // unmask public value
    }

    for(int j=0; j<32; j++) { // batch of 32 r
      if( ((valid[0] >>j)&1) == 0 )
        continue; // r_j is >= N -> next r
    
      for(int s=0; s<NSHARES; s++) {
        for(int k=0; k<LOG_N; k++) {
          r32[s][k] = (rand[s][r+k] >>j) &1; // copy r_j (bitsliced domain)
          r32[s][k] = ~(r32[s][k] + 0xffffffff); // extend single bit to 32bit
        }
      }

      for(c=0; c<i; c+=32) {
        masked_bs_eq(res, 1, &r32[0][0], LOG_N, &p[0][c], W_PADDED, LOG_N); // check for collision
        for(int s=1; s<NSHARES; s++) {
          res[0] ^= res[s]; // unmask public value
        }
        
        index = i &31; // mod 32 -> bit within 32bit register where r gets stored bitsliced

        if(c+31>i-1) // we included comparisons with the initialized 0 above i-1
          res[0] &= 0xffffffff >> (32-index); // mask out these comparisons

        if (res[0] != 0)
          break;  // collision found
      }
      if(res[0] != 0)
        continue; // collision found -> next r

      // no collision:
      if(c > i) 
        c-=32; // c now points to nearest multiple of 32 <= i;

      for(int s=0; s<NSHARES; s++) {
        for(int b=0; b<LOG_N; b++) {
          p[s][c+b] &=  0xffffffff ^ (1<<index); //set bit to zero
          p[s][c+b] |=  (r32[s][b]&1)<<index; //set bit according to r 
        }
      }
  
      i++;
      if(i >= W)
        break;  // stop 32er batch
    }
    r += 32;
  }
  transpose32_array(&p[0][0], W_PADDED); // unbitslice
  return 0;
}
