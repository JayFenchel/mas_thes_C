#include "include/hhmpcalg.h"


void form_Y(real_t Y[],
            const real_t Phi[],
            const uint32_t T,
            const real_t A[], const uint32_t n,
            const real_t B[], const uint32_t m)
{
    
}

real_t *form_Y11(const real_t B[], const uint32_t n, const uint32_t m,
                 const real_t R0_I[],
                 const real_t Q1[])
{
    real_t BT[n*m]; /*TODO Speicher allozieren?*/
    real_t hilf1[m]; /*TODO Speicher allozieren?*/
    real_t hilf2[n*m]; /*TODO Speicher allozieren?*/
    real_t sol[n*n];
    uint32_t i;
    mpcinc_mtx_transpose(BT, B, n, m);
    for (i = 0; i < n; i++)
    {
        fwd_subst(hilf1, R0_I, m, BT+i*m);
        bwd_subst(hilf2+i*m, R0_I, m, hilf1);
    }
    mpcinc_mtx_multiply_mtx_mtx();
    mpcinc_mtx_add();
}

void setBlock(real_t mtx[], const uint32_t dim,
              const real_t blk[], const uint32_t s,
              const uint32_t fst)
{
    uint32_t i, j; /* loop counters */
    
    for (i = 0; i < s; i++){
        for (j = 0; j < s; j++){
            mtx[(i+fst)*dim+j+fst] = blk[i*s+j];
        }
    }
}

void getBlock(real_t blk[],
              const real_t mtx[], const uint32_t dim,
              const uint32_t fst, const uint32_t s)
{
    uint32_t i, j; /* loop counters */
    
    for (i = 0; i < s; i++){
        for (j = 0; j < s; j++){
            blk[i*s+j] = mtx[(i+fst)*dim+j+fst];
        }
    }
}