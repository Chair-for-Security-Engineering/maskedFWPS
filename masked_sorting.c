#include <stdint.h>
#include <stddef.h>
#include "gadgets.h"
#include "params.h"
#include "masked_fwsampling.h"
#include "masked_representations.h"

#define SORTN PAD64(N)

/*
 * a[0] holds the low bits of the first share
 * ...
 * a[31] holds the MSB of the first share
 * a[SORTN] holds the LSBs of the second share
 * ...
 * a[SORTN+31] holds the MSBs of the second share
 * and so on
 */
static void bscmp(uint32_t *res, const uint32_t *a, const uint32_t stride_a, const uint32_t *b, const uint32_t stride_b)
{
  size_t n,i;
  uint32_t xor[NSHARES], tmp[NSHARES];
  
  for (n = 0; n < NSHARES; n++)
  {
    res[n] = 0;
  }
  
  for (i = 2; i < 32; i++) // lowest 2bit do not belong to sorted value
  {
    masked_xor(NSHARES, xor, 1, a + i, stride_a, b + i, stride_b);
    masked_xor(NSHARES, tmp, 1, res, 1, b + i, stride_b);
    masked_and(NSHARES, tmp, 1, tmp, 1, xor, 1);
    masked_xor(NSHARES, res, 1, res, 1, tmp, 1);
  }
}

static void bsint32_MINMAX(uint32_t *a, uint32_t *b, uint32_t direction, uint32_t stride_b)
{
  size_t n,i;
  uint32_t res[NSHARES], xor[NSHARES], mux[NSHARES][32];
  // res = 1 iff b > a
  bscmp(res, a, SORTN, b, stride_b);

  //flip res based on sorting direction
  res[0] = res[0] ^ direction;
  res[0] = ~res[0]; //flip res for faster cmov

  // now swap conditionally based on res
  for (i = 0; i < 32; i++)
  {
    masked_xor(NSHARES, xor, 1, &a[i], SORTN, &b[i], stride_b);
    masked_and(NSHARES, &mux[0][i], 32, xor, 1, res, 1);
  }

  // now, mux will contain a if b > a, else b
  for (n = 0; n < NSHARES; n++)
  {
    for (i = 0; i < 32; i++)
    {
      a[n * SORTN+i] ^= mux[n][i];
      b[n * stride_b+i] ^= mux[n][i];
    }
  }
}

static int isPowerOfTwo (unsigned int x)
{
  return ((x != 0) && !(x & (x - 1)));
}

static int greatestPowerOfTwoLessThan(int n)
{
    int k=1;
    while (k<n)
        k=k<<1;
    return k>>1;
}

static int greatestPowerOfTwoLessThanEqual(int n)
{
    if(isPowerOfTwo(n))
        return n;
    else
      return greatestPowerOfTwoLessThan(n);
}

static uint32_t rotr32a(uint32_t x, uint32_t n)
{
  return (x<<(32-n)) | (x>>n) ;
}

static uint32_t rotl32a(uint32_t x, uint32_t n)
{
  return (x>>(32-n)) | (x<<n) ;
}

// converts two consecutive 32bit words for a minmax step in bitonic sort with p<32
// so that 32 minmax operations can be done in parallel via bitslicing
static void bsconvertsub32(uint32_t *x, int current_distance_log, int target_distance_log) {
  uint32_t t;

  const uint32_t mask[5] = {0x55555555, 0x33333333, 
                0x0F0F0F0F, 0x00FF00FF, 0x0000FFFF};

  for(int c=0; c<32; c++) {
    for(int i=current_distance_log-1; i >= target_distance_log; i--) {
      t = rotr32a(x[c] &~mask[i], 1<<i);
      x[c] = (x[c] & mask[i]) | rotl32a(x[c+32] & mask[i], 1<<i);
      x[c+32] = t | (x[c+32] & ~mask[i]);    
    }
  }
}

static void bsbackconvertsub32(uint32_t *x, int current_distance_log, int target_distance_log) {
   uint32_t t;

  const uint32_t mask[5] = {0x55555555, 0x33333333, 
                0x0F0F0F0F, 0x00FF00FF, 0x0000FFFF};
  
  for(int c=0; c<32; c++) {
    for(int i=current_distance_log; i < target_distance_log; i++) {
      t = rotr32a(x[c] &~mask[i], 1<<i);
      x[c] = (x[c] & mask[i]) | rotl32a(x[c+32] & mask[i], 1<<i);
      x[c+32] = t | (x[c+32] & ~mask[i]);    
    }
  }
}

static void bsconvertsub32poly(uint32_t *x, int n, int current_distance, int target_distance) {
  int current_distance_log = 0;
  int target_distance_log = 0;

  if(target_distance>=32 && current_distance >=32)
    return;
  if(target_distance>32)
    target_distance = 32;
  if(current_distance>32)
    current_distance = 32;

  while(1<<current_distance_log < current_distance)
    current_distance_log++;
  while(1<<target_distance_log < target_distance)
    target_distance_log++;

  if (target_distance < current_distance) {
    for (int s = 0; s < NSHARES; s++) {
      for (int i = 0; i < n - 63; i+=64) {
        bsconvertsub32(&x[s*SORTN + i], current_distance_log, target_distance_log);
      }
    } 
  }
  else if(target_distance > current_distance) {
    for (int s = 0; s < NSHARES; s++) {
      for (int i = 0; i < n - 63; i+=64) {
        bsbackconvertsub32(&x[s*SORTN + i], current_distance_log, target_distance_log);
      }
    }  
  }
}

//expects already bitsliced polynomial
static void masked_bitonic_merge(uint32_t *x, int n, int k, int direction)
{
  uint32_t directionmask;
  int i,j,old_j=32;

  for (j=k; j>0; j=j>>1) {
    bsconvertsub32poly(x, n, old_j, j); // from/to minmax within one register

    if(j>=32) { // pairs do not share registers
      for (i=0; i+j+31<n; i+=32) {
        int ij=i^j;
        if ((ij)>i) {
          directionmask = !(i&k<<1) ? 0x00000000 : 0xFFFFFFFF;
          if(direction)
            directionmask = ~directionmask;
          bsint32_MINMAX(&x[i], &x[ij], directionmask, SORTN);
        }
      }
    }
    else if (j<32) { // pairs do share registers -> seperate them with bsconvertsub32poly function -> thus we need to operate on 64 consecutive elements simulatenously  
      switch (k) {
        case 16: directionmask = 0xFFFF0000; break;
        case  8: directionmask = 0xFF00FF00; break;
        case  4: directionmask = 0xF0F0F0F0; break;
        case  2: directionmask = 0xCCCCCCCC; break;
        case  1: directionmask = 0xAAAAAAAA; break;
        default: directionmask = 0x00000000; break;
      } 
    
      for (i=0; i+63 < n; i+= 64) {
        if(k > 16) {
            directionmask = !(i&k<<1) ? 0x00000000 : 0xFFFFFFFF;
        }
        if(direction)
          directionmask = ~directionmask;  
        bsint32_MINMAX(&x[i], &x[i+32], directionmask, SORTN);
      }
    }
    old_j = j;
  }
  bsconvertsub32poly(x, n, old_j, 32); // from/to minmax within one register
}

//expects already bitsliced polynomial
static void masked_bitonic_sort_pow2(uint32_t *x, int n, int direction)
{
  for (int k=2; k<=n; k=2*k) {
    masked_bitonic_merge(x, n, k>>1, direction);
  }
}

// Very efficient if n is power of two.
// Not the most efficient approach for bitonic sort with an n that is not a power of two,
// we cant reach nlogn because n=64 is the lowest we can get down to while keeping the code somewhat simple.
// Slightly more efficient variant would be possible though.
void masked_bitonic_sort(uint32_t x[NSHARES][SORTN])
{
  transpose32_array(&x[0][0], SORTN); // bitslicing

  int pows2[10] = {0}; //hardcoded indirect limit for size of N
  int n = SORTN;
  int m = 0;
  int i = 0;

  while(n>0) {
    m = greatestPowerOfTwoLessThanEqual(n);
    pows2[i] = m;
    i++;
    n = n-m;
  }

  int32_t *b = &x[0][0] +SORTN-pows2[--i];
  masked_bitonic_sort_pow2(b, pows2[i], 0);
  int n_merge = pows2[i];

  while(i>0) {
    b = b-pows2[--i];
    masked_bitonic_sort_pow2(b, pows2[i], 1);
    n_merge += pows2[i];
    masked_bitonic_merge(b, n_merge, greatestPowerOfTwoLessThan(n_merge), 0);
  }
  
  transpose32_array(&x[0][0], SORTN); // unbitscliing
}

