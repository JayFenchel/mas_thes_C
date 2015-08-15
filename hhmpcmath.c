#include "include/hhmpcmath.h"
#include <hhmpctestfunc.h>


void fwd_subst(real_t sol[],
               const real_t mtx[], const uint32_t dim,
               const real_t vec[])
{
    uint32_t i, j; /* loop counters */
    
    printf("%f, %f, %f, %f\n", mtx[0], mtx[1], mtx[2], mtx[3]);
    
    for (i = 0; i < dim; i++){
        sol[i] = vec[i];
        for (j = 0; j < i; j++){
            sol[i] -= mtx[i*dim+j]*sol[j];
        }
        printf("%f %d\n",mtx[1], 1);
        sol[i] /= mtx[i*dim+j];
    }
}

void bwd_subst(real_t sol[],
               const real_t mtx[], const uint32_t dim,
               const real_t vec[])
{
    uint32_t i, j; /* loop counters */
    
    for (i = 0; i < dim; i++){
        sol[dim-i-1] = vec[dim-i-1];
        for (j = 0; j < i; j++){
            sol[dim-i-1] -= mtx[dim*dim-i*dim-j-1]*sol[dim-j-1];
        }
        printf("%f %d\n",mtx[1], 1);
        sol[dim-i-1] /= mtx[dim*dim-i*dim-j-1];
    }
}

uint32_t mtx_cmp(const real_t mtxA[], const real_t mtxB[], real_t dim)
{
    uint32_t i, cmp = 0;
    for (i = 0; i<dim; i++){
        cmp += (mtxA[i] != mtxB[i]) ? 1 : 0;
    }
    return cmp;
}

real_t mtx_out(real_t mtx[], uint32_t rows, uint32_t cols, real_t vec[]){
    real_t out[2];
    mpcinc_mtx_multiply_mtx_vec(out, mtx, vec, rows, cols);
    return out[0];
}

int simple_sum(int a, int b) {

    return a + b;
}