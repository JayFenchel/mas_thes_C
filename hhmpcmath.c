#include "include/hhmpcmath.h"


void cholesky(real_t sol[],
              const real_t mtx[], const uint32_t dim)
{
    uint32_t i, j, k; /* loop counters */
    
    for (i = 0; i < dim*dim; i++){
        sol[i] = mtx[i];
    }
    
    for (i = 0; i < dim; i++){
//         sol[i*dim+i]+=0.0001;  /* Immernoch benötigt, wenn nicht condensierte Probleme gelöst werden */
        for (j = 0; j < i; j++){
            sol[i*dim+i] -= sol[i*dim+j]*sol[i*dim+j];
        }
//         printf("HIER sqrt\n");
        if (sol[i*dim+i]<0) printf("WARNUNG! %f\n",sol[i*dim+i] ); 
        sol[i*dim+i] = sqrtf(smpl_abs(sol[i*dim+i]));  /* I think it's better if algo fails! */
//         sol[i*dim+i] = smpl_sqrt(smpl_abs(sol[i*dim+i]), 2.);
        
//         printf("HIER after sqrt\n");
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
    uint32_t i, cmp = -1;
    for (i = 0; i<dim; i++){
        cmp = smpl_abs((mtxA[i] - mtxB[i])) <= acc ? cmp : i;
    }
    return cmp;
}

real_t smpl_sqrt(real_t r, real_t e)
{
    if (r != r) {; /*printf("r < 0 :%f\n", r);*/}
    if  (smpl_abs(r/e - e) < SQRT_ACC)
        return e;
    else
        return smpl_sqrt(r, (e + r/e) / 2);
}

real_t nth_root(real_t A, uint32_t n){
    if (n == 0) return 1.;
    const uint32_t K = 3;
    uint32_t k;
//     printf("n=%d\n",n);
    real_t x[] = {1., 0., 0.};
    for (k = 0; k < K - 1; k++){
//         printf("in sp3\n");
        x[k + 1] = (1./n) * ((n - 1)*x[k] + A/smpl_pow(&x[k], n - 1));
    }
//     printf("f=%f\n",x[K-1]);
    return x[K-1];
}

real_t smpl_pow(real_t *b, real_t e)
{
//     real_t half_pow, e_halbe;
//     printf("in sp1\n");
//     printf("in sp2 %f\n", e);
//     if (e != e){
//          printf("in sp6 %f\n", e);
//         return e;}
//     if (e == 0.){
//         return 1.;
//         
//     }else if (e < 0){
//         printf("in sp5 %f\n", e_halbe);
//         return 1 / smpl_pow(b, -e);
//     }else if (e > 0. && e < 1){
// //         printf("%f\n", (1./((uint32_t)(1/(1./((uint32_t)(1/e)) - e))) - (1./((uint32_t)(1/e)) - e)));
//         printf("in sp\n");
//         return nth_root(b[0], 1./e) / ( nth_root(b[0], 1./(1./((uint32_t)(1/e)) - e)) / nth_root(b[0], 1./(1./((uint32_t)(1/(1./((uint32_t)(1/e)) - e))) - (1./((uint32_t)(1/e)) - e))) );
//     }else if ((uint32_t)e % 2 == 0){
//         printf("in sp2\n");
//         e_halbe = e/2;
//         
//         half_pow = smpl_pow(b, e_halbe);
//         
//         return half_pow * half_pow;
//     }else{
//         return b[0] * smpl_pow(b, e - 1);
//     }
    return expf(e);
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