#ifndef MASKED_FWSAMPLING__H
#define MASKED_FWSAMPLING__H

#include "params.h"
#include "masked_representations.h"

// BIN2TRI
void masked_binary_to_trinary(uint32_t signs[NSHARES][(N+31)/32], const uint32_t poly[NSHARES][(N+31)/32]);
void masked_binary_to_trinary_fixed(uint32_t signs[NSHARES][(N+31)/32], const uint32_t poly[NSHARES][(N+31)/32], const uint8_t target_weight);

// I2C
void bs_coeff_to_regular(uint32_t poly[NSHARES][N_PADDED/32], uint32_t index[NSHARES][W_PADDED]);
// expects non bitsliced inputs
void bs_coeff_to_regular_fast(uint32_t poly[NSHARES][N_PADDED/32], uint32_t index[NSHARES][W_PADDED]);

// CMP
uint8_t sample_fixed(uint32_t out[NSHARES][(N+31)/32]);

// FISHER-YATES
//expects randomness in p for first W entries
void bs_fisher_yates(uint32_t p[NSHARES][W_PADDED]);

// REJECT
//expects non bitsliced input
//expects rand filled with randomness
int bs_ct_reject(uint32_t p[NSHARES][W_PADDED], uint32_t rand[NSHARES][BOUND]);
int bs_reject(uint32_t p[NSHARES][W_PADDED], uint32_t rand[NSHARES][B_PADDED]);

// REPAND
void repand_sample_fixed(uint32_t A[NSHARES][(N+31)/32], const uint16_t target_weight); // we assume A to be initialized to zero!

#endif
