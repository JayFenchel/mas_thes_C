#include "include/hhmpcmath.h"


void cholesky(real_t sol[],
              const real_t mtx[], const uint32_t dim)
{
    uint32_t i, j, k; /* loop counters */
    
    for (i = 0; i < dim*dim; i++){
        sol[i] = mtx[i];
    }
    
    for (i = 0; i < dim; i++){
        for (j = 0; j < i; j++){
            sol[i*dim+i] -= sol[i*dim+j]*sol[i*dim+j];
        }
        sol[i*dim+i] = smpl_sqrt(smpl_abs(sol[i*dim+i]), 2.);
        
        for (j = i+1; j < dim; j++){
            for (k = 0; k < i; k++){
                sol[j*dim+i] -= sol[j*dim+k]*sol[i*dim+k];
            }
            sol[j*dim+i] /= sol[i*dim+i];             
        }        
    }
    for (i = 0; i < dim; i++){
        for (j = i+1; j < dim; j++){
            sol[i*dim+j] = 0.;
        }
    }
}

void fwd_subst(real_t sol[],
               const real_t mtxA[], const uint32_t dim,
               const real_t mtxB[], const uint32_t colsB)
{
    uint32_t i, j, k; /* loop counters */
    for (k = 0; k < colsB; k++){
        for (i = 0; i < dim; i++){
            sol[colsB*i+k] = mtxB[colsB*i+k];
            for (j = 0; j < i; j++){
                sol[colsB*i+k] -= mtxA[i*dim+j]*sol[colsB*j+k];
            }
            sol[colsB*i+k] /= mtxA[i*dim+j];
        }
    }
    
}

void bwd_subst(real_t sol[],
               const real_t mtxA[], const uint32_t dim,
               const real_t mtxB[], const uint32_t colsB)
{
    uint32_t i, j, k; /* loop counters */
    for (k = colsB; k > 0; k--){
        for (i = 0; i < dim; i++){
            sol[colsB*(dim-i)-k] = mtxB[colsB*(dim-i)-k];
            for (j = 0; j < i; j++){
                sol[colsB*(dim-i)-k] -= mtxA[dim*dim-i*dim-j-1]*sol[colsB*(dim-j)-k];
            }
            sol[colsB*(dim-i)-k] /= mtxA[dim*dim-i*dim-j-1];
        }
    }
}

uint32_t mtx_cmp(const real_t mtxA[], const real_t mtxB[], real_t dim, real_t acc)
{
    uint32_t i, cmp = 0;
    for (i = 0; i<dim; i++){
        cmp += smpl_abs((mtxA[i] - mtxB[i])) > acc ? 1 : 0;
        if (smpl_abs((mtxA[i] - mtxB[i])) > acc)
        printf("%d\n", i);
    }
    return cmp;
}

real_t smpl_sqrt(real_t r, real_t e)
{
    if  (smpl_abs(r/e - e) < SQRT_ACC)
        return e;
    else
        return smpl_sqrt(r, (e + r/e) / 2);
}

real_t smpl_abs(real_t x)
{
    return (x < 0) ? -x : x;
}
    

real_t mtx_out(real_t mtx[], uint32_t rows, uint32_t cols, real_t vec[]){
    real_t out[2];
    mpcinc_mtx_multiply_mtx_vec(out, mtx, vec, rows, cols);
    return out[0];
}

int simple_sum(int a, int b) {

    return a + b;
}