#include "include/hhmpcalg.h"


void solve_sysofleq(real_t delta_z[], real_t delta_v[],
                    const real_t Phi[],
                    const real_t rd[], const real_t rp[],
                    const real_t C[], const real_t A[], const real_t B[],
                    const uint32_t n, const uint32_t m, const uint32_t T)
{
    real_t L_Phi_blocks[m*m + (T-1)*(n+m)*(n+m) + n*n]; /*blocks discribed in paper*/
    real_t L_Phi[T*(n+m)*T*(n+m)];
    real_t beta[T*n];
    real_t C_T[T*(n+m) * T*n];
    real_t Y[T*n*T*n];
    zeroes(Y, T*n*T*n);
    zeroes(L_Phi, T*(n+m)*T*(n+m));


    form_Y(Y, L_Phi_blocks, Phi, T, A, n, B, m);

    /*Ohne Schleife klappt es so nur für T = 3*/
    setBlock(L_Phi, T*(n+m), L_Phi_blocks, m, m, 0, 0);
    setBlock(L_Phi, T*(n+m), L_Phi_blocks+m*m, n+m, n+m, m, m);
    setBlock(L_Phi, T*(n+m), L_Phi_blocks+m*m+1*(n+m)*(n+m), n+m, n+m, m+1*(n+m), m+1*(n+m));
    setBlock(L_Phi, T*(n+m), L_Phi_blocks+m*m+2*(n+m)*(n+m), n, n, m+2*(n+m), m+2*(n+m));
    
    mpcinc_mtx_transpose(C_T, C, T*n, T*(n+m));
    
    form_beta(beta, L_Phi_blocks, rd, rp, T, C, n, m);
    form_delta_v(delta_v, Y, beta, T, n);
    form_delta_z(delta_z, delta_v, L_Phi, rd, C_T, T, n, m);   
}

void form_delta_z(real_t delta_z[],
                  const real_t delta_v[],
                  const real_t L_Phi[],
                  const real_t rd[],
                  const real_t C_T[],
                  const uint32_t T, const uint32_t n, const uint32_t m)
{
    real_t help1[T*(n+m)];
    real_t help2[T*(n+m)];
    real_t L_Phi_T[T*(n+m)*T*(n+m)];
    uint32_t i;
    mpcinc_mtx_multiply_mtx_vec(help1, C_T, delta_v, T*(n+m), T*n);
    for (i = 0; i < T*(n+m); i++)
        help1[i] = -help1[i];
    mpcinc_mtx_substract(help2, help1, rd, T*(n+m), 1);
    mpcinc_mtx_transpose(L_Phi_T, L_Phi, T*(n+m), T*(n+m));
    fwd_subst(help1, L_Phi, T*(n+m), help2, 1);
    bwd_subst(delta_z, L_Phi_T, T*(n+m), help1, 1);
    
    
}


void form_delta_v(real_t delta_v[],
                  const real_t Y[],
                  const real_t beta[],
                  const uint32_t T, const uint32_t n)
{
    real_t L_Y[T*n*T*n];
    real_t L_Y_T[T*n*T*n];
    real_t help1[T*n];
    real_t mbeta[T*n];
    uint32_t i;
    
    for (i = 0; i < T*n; i++)
        mbeta[i] = -beta[i];
    
    
    cholesky(L_Y, Y, T*n);
    mpcinc_mtx_transpose(L_Y_T, L_Y, T*n, T*n);
    fwd_subst(help1, L_Y, T*n, mbeta, 1);
    bwd_subst(delta_v, L_Y_T, T*n, help1, 1);
    
}

void form_beta(real_t beta[],
               const real_t L_Phi_blocks[],
               const real_t rd[], const real_t rp[],
               const uint32_t T,
               const real_t C[], const uint32_t n, const uint32_t m
               /*const real_t A[], const uint32_t n,
               const real_t B[], const uint32_t m*/)
{
    real_t L_Phi[T*(n+m)*T*(n+m)];
    real_t L_Phi_T[T*(n+m)*T*(n+m)];
    real_t help1[T*(n+m)];
    real_t help2[T*(n+m)];
    real_t help3[T*n];
    
    uint32_t g;
    for(g = 0; g < T*(n+m)*T*(n+m); g++)
        L_Phi[g] = 0;
    
    /*Ohne Schleife klappt es so nur für T = 3*/
    setBlock(L_Phi, T*(n+m), L_Phi_blocks, m, m, 0, 0);
    setBlock(L_Phi, T*(n+m), L_Phi_blocks+m*m, n+m, n+m, m, m);
    setBlock(L_Phi, T*(n+m), L_Phi_blocks+m*m+1*(n+m)*(n+m), n+m, n+m, m+1*(n+m), m+1*(n+m));
    setBlock(L_Phi, T*(n+m), L_Phi_blocks+m*m+2*(n+m)*(n+m), n, n, m+2*(n+m), m+2*(n+m));
    
    
    mpcinc_mtx_transpose(L_Phi_T, L_Phi, T*(n+m), T*(n+m));
    fwd_subst(help1, L_Phi, T*(n+m), rd, 1);
    bwd_subst(help2, L_Phi_T, T*(n+m), help1, 1);
    mpcinc_mtx_multiply_mtx_vec(help3, C, help2, T*n, T*(n+m));
    mpcinc_mtx_substract(beta, help3, rp, T*n, 1);
}

void form_Y(real_t Y[], real_t L_Phi[],
            const real_t Phi[],
            const uint32_t T,
            const real_t A[], const uint32_t n,
            const real_t B[], const uint32_t m)
{
    uint32_t i, j;
    real_t R0_I[m*m];
    real_t PhiBlock[(n+m)*(n+m)];
    real_t PhiBlock_C_T[(n+m)*(n+m)];
    real_t hilf1[(n+m)*(n+m)];
    real_t hilf2[(n+m)*(n+m)];
    real_t last_hilf2[(n+m)*(n+m)];
    real_t Q1[n*n];
    real_t Y11[n*n];
    real_t eye1[(n+m)*(n+m)];
    real_t eye2[n*n];
    eye(eye1, n+m);
    eye(eye2, n);
    
    real_t A_T[n*n];
    real_t B_T[m*n];
    
    real_t Q_I[n*n];
    real_t Q_I_C_T[n*n];
    real_t hilf4[n*n];
    
    mpcinc_mtx_transpose(B_T, B, n, m);
    mpcinc_mtx_transpose(A_T, A, n, n);
        
    real_t A_T_B_T[n*(n+m)];
    for (i = 0; i < n*n; i++)
        A_T_B_T[i] = A_T[i];
    for (i = 0; i < n*m; i++)
        A_T_B_T[n*n+i] = B_T[i];

    real_t A_B[(n+m)*n];
    mpcinc_mtx_transpose(A_B, A_T_B_T, n+m, n);
    
    for (i = 0; i < T; i++){
    if (i == 0){
    
    
    getBlock(R0_I, Phi, T*(n+m), 0, 0, m, m);
    cholesky(L_Phi, R0_I, m);
    getBlock(PhiBlock, Phi, T*(n+m), m+i*(n+m), m+i*(n+m), n+m, n+m);
    cholesky(L_Phi+m*m, PhiBlock, n+m);
    fwd_subst(hilf1, L_Phi+m*m, n+m, eye1, n+m);
    mpcinc_mtx_transpose(PhiBlock_C_T, L_Phi+m*m, n+m, n+m);
    bwd_subst(hilf2, PhiBlock_C_T, n+m, hilf1, n+m);
    getBlock(Q1, hilf2, n+m, 0, 0, n, n);
    form_Y11(Y11, B, n, m, L_Phi, Q1);
    setBlock(Y, T*n, Y11, n, n, 0, 0);
    
    }
    
    if (i > 0){
        for (j = 0; j < (n+m)*(n+m); j++)
            last_hilf2[j] = hilf2[j];
        if (i == T-1){
            getBlock(Q_I, Phi, T*(n+m), m+i*(n+m), m+i*(n+m), n, n);
            cholesky(L_Phi+m*m+i*(n+m)*(n+m), Q_I, n);
            fwd_subst(hilf4, L_Phi+m*m+i*(n+m)*(n+m), n, eye2, n);
            mpcinc_mtx_transpose(Q_I_C_T, L_Phi+m*m+i*(n+m)*(n+m), n, n);
            bwd_subst(Q1, Q_I_C_T, n, hilf4, n);
            
            
        }else{
        getBlock(PhiBlock, Phi, T*(n+m), m+i*(n+m), m+i*(n+m), n+m, n+m);
        cholesky(L_Phi+m*m+i*(n+m)*(n+m), PhiBlock, n+m);
        fwd_subst(hilf1, L_Phi+m*m+i*(n+m)*(n+m), n+m, eye1, n+m);
        mpcinc_mtx_transpose(PhiBlock_C_T, L_Phi+m*m+i*(n+m)*(n+m), n+m, n+m);
        bwd_subst(hilf2, PhiBlock_C_T, n+m, hilf1, n+m);
        getBlock(Q1, hilf2, n+m, 0, 0, n, n);
            
        }
        
        form_Yii(Y11, A_B, n, n+m, last_hilf2, n+m, n+m, A_T_B_T, n+m, n, Q1);
        setBlock(Y, T*n, Y11, n, n, i*n, i*n);
        
    }
    
    if (i < T-1){
        form_Y_i_ip1(Y11, A_T, n, B_T, m, hilf2);
        setBlock(Y, T*n, Y11, n, n, i*n, (i+1)*n);
        real_t Y11_T[n*n];
        mpcinc_mtx_transpose(Y11_T, Y11, n, n);
        setBlock(Y, T*n, Y11_T, n, n, (i+1)*n, i*n);
    }   
    
    }
}

void form_Yii(real_t solution[],
              const real_t A_B[], const uint32_t rowsA, const uint32_t colsA,
              const real_t hilf2[], const uint32_t rowsB, const uint32_t colsB,
              const real_t A_T_B_T[], const uint32_t rowsC, const uint32_t colsC,
              const real_t Qi[])
{
    real_t help1[rowsA*(colsA+colsB)];
    real_t help2[colsA*rowsA];
    
    mpcinc_mtx_multiply_mtx_mtx(help1, hilf2, A_T_B_T, rowsB, colsB, colsC);
    mpcinc_mtx_multiply_mtx_mtx(help2, A_B, help1, rowsA, colsA, colsC);
    mpcinc_mtx_add(solution, help2, Qi, rowsA, rowsA);
}

void form_Y_i_ip1(real_t sol[],
                  const real_t A_T[], const uint32_t n,
                  const real_t B_T[], const uint32_t m,
                  const real_t hilf2[])
{
    uint32_t i;
    real_t mA_T_B_T[n*(n+m)];
    for (i = 0; i < n*n; i++)
        mA_T_B_T[i] = -A_T[i];
    for (i = 0; i < n*m; i++)
        mA_T_B_T[n*n+i] = -B_T[i];
    
    mpcinc_mtx_multiply_mtx_mtx(sol, hilf2, mA_T_B_T, n+m, n+m, n);

}

void form_Y11(real_t sol[],
              const real_t B[], const uint32_t n, const uint32_t m,
              const real_t R0_C[],
              const real_t Q1[])
{
    real_t BT[m*n]; /*TODO Speicher allozieren?*/
    real_t hilf1[m*n]; /*TODO Speicher allozieren?*/
    real_t hilf2[m*n]; /*TODO Speicher allozieren?*/
    real_t R0_C_T[m*m];
    real_t hilf3[n*n]; /*TODO Speicher allozieren?*/
    mpcinc_mtx_transpose(BT, B, n, m);
    fwd_subst(hilf1, R0_C, m, BT, n);
    mpcinc_mtx_transpose(R0_C_T, R0_C, m, m);
    bwd_subst(hilf2, R0_C_T, m, hilf1, n);
    mpcinc_mtx_multiply_mtx_mtx(hilf3, B, hilf2, n, m, n);
    mpcinc_mtx_add(sol, hilf3, Q1, n, n);
}

void setBlock(real_t mtx[], const uint32_t dim,
              const real_t blk[], const uint32_t s_row, const uint32_t s_col,
              const uint32_t row_fst, const uint32_t col_fst)
{
    uint32_t i, j; /* loop counters */
    
    for (i = 0; i < s_col; i++){
        for (j = 0; j < s_row; j++){
            mtx[(i+row_fst)*dim+j+col_fst] = blk[i*s_row+j];
        }
    }
}

void getBlock(real_t blk[],
              const real_t mtx[], const uint32_t dim,
              const uint32_t row_fst, const uint32_t col_fst,
              const uint32_t s_row, const uint32_t s_col)
{
    uint32_t i, j; /* loop counters */
    
    for (i = 0; i < s_col; i++){
        for (j = 0; j < s_row; j++){
            blk[i*s_row+j] = mtx[(i+row_fst)*dim+j+col_fst];
        }
    }
}