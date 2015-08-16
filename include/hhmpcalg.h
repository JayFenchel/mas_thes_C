#ifndef HHMPCALG_H
#define HHMPCALG_H

#include "arithmetic.h"
#include "hhmpcmath.h"
#include "mpcincmtxops.h"


void form_Y(real_t mtxY[],
            const real_t mtxPhi[],
            const uint32_t horizon,
            const real_t mtxA[], const uint32_t dimA,
            const real_t mtxB[], const uint32_t colsB);

real_t *form_Y11(const real_t B[], const uint32_t rowsB, const uint32_t colsB,
                 const real_t R0[],
                 const real_t Q1[]);

void setBlock(real_t mtx[], const uint32_t dimension, 
              const real_t block[], const uint32_t size,
              const uint32_t first);

void getBlock(real_t block[], 
              const real_t mtx[], const uint32_t dimension,
              const uint32_t first,const uint32_t size);


#endif