#ifndef HHMPCALG_H
#define HHMPCALG_H

#include "arithmetic.h"
#include "hhmpcmath.h"
#include "mpcincmtxops.h"
#include "hhmpcusefull.h"


void form_delta_z(real_t delta_z[],
                  const real_t delta_v[],
                  const real_t mtxL_Phi[],
                  const real_t rd[],
                  const real_t C[],
                  const uint32_t T, const uint32_t n, const uint32_t m);

void form_delta_v(real_t delta_v[],
                  const real_t Y[], 
                  const real_t beta[],
                  const uint32_t T, const uint32_t n);

void form_beta(real_t beta[],
               const real_t mtxL_Phi[],
               const real_t rd[], const real_t rp[],
               const uint32_t horizon,
               const real_t mtxA[], const uint32_t dimA,
               /*const real_t mtxB[],*/ const uint32_t colsB);

/*returns mtxL_Phi the cholesky factorization of Phi
 */
void form_Y(real_t mtxY[], real_t mtxL_Phi[],
            const real_t mtxPhi[],
            const uint32_t horizon,
            const real_t mtxA[], const uint32_t dimA,
            const real_t mtxB[], const uint32_t colsB);

void form_Yii(real_t solution[],
              const real_t A[], const uint32_t rowsA, const uint32_t colsA,
              const real_t B[], const uint32_t rowsB, const uint32_t colsB,
              const real_t C[], const uint32_t rowsC, const uint32_t colsC,
              const real_t Q[]);

void form_Y_i_ip1(real_t solution[],
                  const real_t A[], const uint32_t rowsA, const uint32_t colsA,
                  const real_t B[], const uint32_t rowsB, const uint32_t colsB,
                  const real_t Qi_C[]);

void form_Y11(real_t solution[],
              const real_t B[], const uint32_t rowsB, const uint32_t colsB,
              const real_t R0_Cholesky[],
              const real_t Q1[]);

void setBlock(real_t mtx[], const uint32_t dimension, 
              const real_t block[], const uint32_t size_row, const uint32_t size_col,
              const uint32_t first_row, const uint32_t first_col);

void getBlock(real_t block[], 
              const real_t mtx[], const uint32_t dimension,
              const uint32_t first_row, const uint32_t first_col,
              const uint32_t size_row, const uint32_t size_col);


#endif