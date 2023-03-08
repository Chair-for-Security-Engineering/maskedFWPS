#ifndef PARAMS__H
#define PARAMS__H

#define NSHARES 2

#define N 8192
#define W 4
#define BOUND 3

#define LOG_N 13
#define LOG_W 2
#define LOG_BOUND 2

#define PAD32(X) ((((X) + 31)/32)*32)
#define N_PADDED PAD32(N)
#define W_PADDED PAD32(W)
#define B_PADDED PAD32(BOUND)

#define PAD64(X) ((((X) + 63)/64)*64)

#ifndef COEF_NBITS
#define COEF_NBITS 1
#endif

union U64 { 
  uint64_t v64; 
  uint32_t v32[2]; 
};

#endif
