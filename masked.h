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
#ifndef MASKED_H
#define MASKED_H
#include <stdint.h>
#include "randombytes.h"

#ifndef NSHARES
#define NSHARES 2
#endif

#define BSSIZE 32

typedef uint32_t BsBBit[NSHARES];   // dense
typedef uint16_t Coef;
typedef Coef ACoef[NSHARES];            // dense


static inline uint32_t get_random(void) 
{
#if 0 // original ARM code
  while (1) {
    if ((RNG_SR & RNG_SR_DRDY) == 1) { // check if data is ready
      return RNG_DR;
    }
  }
#endif
  uint32_t r;
  randombytes((uint8_t*)&r, sizeof(uint32_t));
  return r;
}

#endif
