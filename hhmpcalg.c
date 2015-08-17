#include "include/hhmpcalg.h"


void form_Y(real_t Y[],
            const real_t Phi[],
            const uint32_t T,
            const real_t A[], const uint32_t n,
            const real_t B[], const uint32_t m)
{
    
}

void form_Y11(real_t sol[],
              real_t B[], const uint32_t n, const uint32_t m,
              const real_t R0_I[],
              const real_t Q1[])
{
    real_t BT[m*n]; /*TODO Speicher allozieren?*/
    real_t hilf1[m*n]; /*TODO Speicher allozieren?*/
    real_t hilf2[m*n]; /*TODO Speicher allozieren?*/
    real_t hilf3[n*n]; /*TODO Speicher allozieren?*/
    mpcinc_mtx_transpose(BT, B, n, m);
    fwd_subst(hilf1, R0_I, m, BT, n);
    bwd_subst(hilf2, R0_I, m, hilf1, n);
    mpcinc_mtx_multiply_mtx_mtx(hilf3, B, hilf2, n, m, n);
    mpcinc_mtx_add(sol, hilf3, Q1, n, n);
}

void setBlock(real_t mtx[], const uint32_t dim,
              const real_t blk[], const uint32_t s,
              const uint32_t row_fst, const uint32_t col_fst)
{
    uint32_t i, j; /* loop counters */
    
    for (i = 0; i < s; i++){
        for (j = 0; j < s; j++){
            mtx[(i+row_fst)*dim+j+col_fst] = blk[i*s+j];
        }
    }
}

void getBlock(real_t blk[],
              const real_t mtx[], const uint32_t dim,
              const uint32_t row_fst, const uint32_t col_fst,
              const uint32_t s)
{
    uint32_t i, j; /* loop counters */
    
    for (i = 0; i < s; i++){
        for (j = 0; j < s; j++){
            blk[i*s+j] = mtx[(i+row_fst)*dim+j+col_fst];
        }
    }
}