#include <stdint.h>
#include <stddef.h>
#include "gadgets.h"
#include "params.h"
#include "masked_fwsampling.h"

// expects non bitsliced inputs
void bs_coeff_to_regular(uint32_t poly[NSHARES][N_PADDED/32], uint32_t index[NSHARES][W_PADDED]) {

  uint32_t p_index[NSHARES][32];

  // init arrays:
  for(int c=0; c<N_PADDED/32; c++) {
    poly[0][c] = 0; // init with 0
  } 

  for(int s=1; s<NSHARES; s++) {
    for(int c=0; c<N_PADDED/32; c++) {
      poly[s][c] = 0; // init with 0
    } 
  }

  for(int s=1; s<NSHARES; s++) {
    for(int c=0; c<32; c++) {
      p_index[s][c] = 0; // public value, no masking
    } 
  }

  uint32_t index_i[NSHARES][32]; // 32 times index i
  uint32_t res[NSHARES]; // result of comparison

  for(int i=0; i<W; i++) { // iterate over all nonzero indices
    // copy index_i 32 times:
    for(int s=0; s<NSHARES; s++) {
      for(int c=0; c<32; c++) {
        index_i[s][c] = index[s][i];
      }
      transpose32(&index_i[s][0]); 
    }
    // iterate over poly, 32 coefficients at a time: 
    for(int c=0; c<N_PADDED; c+=32) {
      for(int j=0; j<32; j++)
        p_index[0][j] = c+j;
      transpose32(&p_index[0][0]);
      // index of coefficient equal to the current nonzero coefficient?
      masked_bs_eq(res, 1, &p_index[0][0], 32, &index_i[0][0], 32, LOG_N);
      // if yes, set coefficient in poly to 1
      masked_bs_sel(res, 1, &poly[0][c/32], N_PADDED/32, res, 1, 1);
    }
  }
  //transpose32_array(&poly[0][0], N_PADDED); // unbitslice result
}

// expects non bitsliced inputs
void bs_coeff_to_regular_fast(uint32_t poly[NSHARES][N_PADDED/32], uint32_t index[NSHARES][W_PADDED]) {

  uint32_t p_index[NSHARES][N_PADDED];

  // init arrays:
  for(int c=0; c<N_PADDED/32; c++) {
    poly[0][c] = 0; // init with 0
  } 

  for(int s=1; s<NSHARES; s++) {
    for(int c=0; c<N_PADDED/32; c++) {
      poly[s][c] = 0; // init with 0
    } 
  }

  for(int c=0; c<N_PADDED; c++) {
    p_index[0][c] = c; // init 
  }
  for(int c=0; c<N_PADDED; c+=32) {
    transpose32(&p_index[0][c]);
  } 

  for(int s=1; s<NSHARES; s++) {
    for(int c=0; c<N_PADDED; c++) {
      p_index[s][c] = 0; // public value, no masking
    } 
  }

  uint32_t index_i[NSHARES][32]; // 32 times index i
  uint32_t res[NSHARES]; // result of comparison

  for(int i=0; i<W; i++) { // iterate over all nonzero indices
    // copy index_i 32 times:
    for(int s=0; s<NSHARES; s++) {
      for(int c=0; c<32; c++) {
        index_i[s][c] = index[s][i];
      }
      transpose32(&index_i[s][0]); 
    }
    // iterate over poly, 32 coefficients at a time: 
    for(int c=0; c<N_PADDED; c+=32) {
      // index of coefficient equal to the current nonzero coefficient?
      masked_bs_eq(res, 1, &p_index[0][c], N_PADDED, &index_i[0][0], 32, LOG_N);
      // if yes, set coefficient in poly to 1
      masked_bs_sel(res, 1, &poly[0][c/32], N_PADDED/32, res, 1, 1);
    }
  }
  //transpose32_array(&poly[0][0], N_PADDED); // unbitslice result
}

