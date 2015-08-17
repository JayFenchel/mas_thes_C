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

void form_Y11(real_t solution[],
              real_t B[], const uint32_t rowsB, const uint32_t colsB,
              const real_t R0_I[],
              const real_t Q1[]);

void setBlock(real_t mtx[], const uint32_t dimension, 
              const real_t block[], const uint32_t size,
              const uint32_t first_row, const uint32_t first_col);

void getBlock(real_t block[], 
              const real_t mtx[], const uint32_t dimension,
              const uint32_t first_row, const uint32_t first_col,
              const uint32_t size);


#endif