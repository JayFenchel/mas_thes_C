#ifndef HHMPCMATH_H
#define HHMPCMATH_H

/* Eulersche Zahl */
#define E 2.71828182845904523536028747135266249775724709369995
#define SQRT_ACC 1e-10

#include "arithmetic.h"
#include "mpcincmtxops.h"


void cholesky(real_t solution[],
              const real_t mtx[], const uint32_t dimension);

/* forward substitution: mtxA * solution = mtxB
 * mtxA has shape (dimA x dimA), lower Dreiecksmatrix
 * sulution and mtxB have shape (dimA x columsB)*/
void fwd_subst(real_t solution[],
               const real_t mtxA[], const uint32_t dimensionA,
               const real_t mtxB[], const uint32_t columsB);

/* backward substitution: mtxA * solution = mtxB
 * mtxA has shape (dimA x dimA), upper Dreiecksmatrix
 * sulution and mtxB have shape (dimA x columsB)*/
void bwd_subst(real_t solution[],
               const real_t mtxA[], const uint32_t dimensionA,
               const real_t mtxB[], const uint32_t columsB);

uint32_t mtx_cmp(const real_t mtxA[], const real_t mtxB[], real_t dim, real_t accuracy);

real_t smpl_sqrt(real_t radikant, real_t sqrt_guess);

real_t nth_root(real_t A, int n);
real_t smpl_pow(real_t base, real_t exp);

real_t smpl_abs(real_t x);

real_t mtx_out(real_t mtx[], uint32_t rows, uint32_t cols, real_t vec[]);

int simple_sum(int a, int b);


#endif