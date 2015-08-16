#ifndef HHMPCMATH_H
#define HHMPCMATH_H

#define SQRT_ACC 1e-10

#include "arithmetic.h"
#include "mpcincmtxops.h"


void cholesky(real_t solution[],
              const real_t mtx[], const uint32_t dimension);

void fwd_subst(real_t solution[],
               const real_t mtx[], const uint32_t dimension,
               const real_t vec[]);

void bwd_subst(real_t solution[],
               const real_t mtx[], const uint32_t dimension,
               const real_t vec[]);

uint32_t mtx_cmp(const real_t mtxA[], const real_t mtxB[], real_t dim);

real_t smpl_sqrt(real_t radikant, real_t sqrt_exp);

real_t smpl_abs(real_t x);

real_t mtx_out(real_t mtx[], uint32_t rows, uint32_t cols, real_t vec[]);

int simple_sum(int a, int b);


#endif