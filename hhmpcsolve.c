#include "include/hhmpcsolve.h"


void solve_sysofleq(real_t delta_z[], real_t delta_v[],
                    const real_t Phi[],
                    const real_t rd[], const real_t rp[],
                    const real_t C[], const real_t *C_T,
                    const real_t A[], const real_t A_T[],
                    const real_t B[], const real_t B_T[],
                    const uint32_t n, const uint32_t m, const uint32_t T,
                    const real_t *eye_nm, const real_t *eye_n,
                    real_t *tmp_optvar_seqlen,
                    real_t *tmp_dual_seqlen,
                    real_t *L_Y, real_t *L_Y_T)
{
    real_t PhiBlock[(n+m)*(n+m)];
    real_t PhiBlock_I[(n+m)*(n+m)];
    real_t PhiBlock_I_last[(n+m)*(n+m)];
    real_t Block_nxn1[n*n];
    real_t Block_nxn2[n*n];
    real_t tmp_optvar_veclenxoptvar_veclen[(n+m)*(n+m)];
    
    real_t L_Phi_blocks[m*m + (T-1)*(n+m)*(n+m) + n*n]; /*blocks discribed in paper*/
    real_t L_Phi_T_blocks[m*m + (T-1)*(n+m)*(n+m) + n*n]; /*blocks discribed in paper*/
    real_t L_Phi[T*(n+m)*T*(n+m)];
    real_t L_Phi_T[T*(n+m)*T*(n+m)];
    real_t beta[T*n];
    real_t Y[T*n*T*n];
    zeroes(Y, T*n*T*n);
    zeroes(L_Phi, T*(n+m)*T*(n+m));
    
    real_t L_Y_blocks[(2*T-1)*n*n];
    real_t L_Y_T_blocks[(2*T-1)*n*n];
    
    real_t A_T_B_T[n*(n+m)];
    uint32_t i;
    for (i = 0; i < n*n; i++)
        A_T_B_T[i] = A_T[i];
    for (i = 0; i < n*m; i++)
        A_T_B_T[n*n+i] = B_T[i];

    real_t A_B[(n+m)*n];
    mpcinc_mtx_transpose(A_B, A_T_B_T, n+m, n);
    
    
    form_Y(L_Y_blocks, L_Y_T_blocks, L_Phi_blocks, L_Phi_T_blocks,
           Phi, T, A_B, A_T_B_T, n, B, B_T, m, eye_nm, eye_n,
           PhiBlock, PhiBlock_I, PhiBlock_I_last,
           Block_nxn1, Block_nxn2, tmp_optvar_veclenxoptvar_veclen);

    /*Ohne Schleife klappt es so nur für T = 3*/
    
    setBlock(L_Phi, T*(n+m), L_Phi_blocks, m, m, 0, 0);
    setBlock(L_Phi, T*(n+m), L_Phi_blocks+m*m, n+m, n+m, m, m);
    setBlock(L_Phi, T*(n+m), L_Phi_blocks+m*m+1*(n+m)*(n+m), n+m, n+m, m+1*(n+m), m+1*(n+m));
    setBlock(L_Phi, T*(n+m), L_Phi_blocks+m*m+2*(n+m)*(n+m), n, n, m+2*(n+m), m+2*(n+m));
    
    mpcinc_mtx_transpose(L_Phi_T, L_Phi, T*(n+m), T*(n+m));
    
    form_beta(beta, L_Phi, L_Phi_T, rd, rp, T, C, n, m);
    form_delta_v(delta_v, tmp_dual_seqlen, L_Y_blocks, L_Y_T_blocks, beta, T, n);
    form_delta_z(delta_z, tmp_optvar_seqlen, delta_v,
                 L_Phi_blocks, L_Phi_T_blocks, rd, C_T, T, n, m);   
}

void form_delta_z(real_t delta_z[],
                  real_t *tmp_optvar_seqlen,
                  const real_t delta_v[],
                  const real_t L_Phi[],
                  const real_t L_Phi_T[],
                  const real_t rd[],
                  const real_t C_T[],
                  const uint32_t T, const uint32_t n, const uint32_t m)
{
    uint32_t i;
    mpcinc_mtx_multiply_mtx_vec(delta_z, C_T, delta_v, T*(n+m), T*n);
    mpcinc_mtx_add_direct(delta_z, rd, T*(n+m), 1);
    mpcinc_mtx_scale_direct(delta_z, -1, T*(n+m), 1);
    
//     fwd_subst(tmp_optvar_seqlen, L_Phi, T*(n+m), delta_z, 1);
//     print_mtx(tmp_optvar_seqlen, T*(n+m), 1);
//     bwd_subst(delta_z, L_Phi_T, T*(n+m), tmp_optvar_seqlen, 1);
    
    fwd_subst(tmp_optvar_seqlen, L_Phi, m, delta_z, 1);
    for (i = 0; i < T-1; i++){
        fwd_subst(tmp_optvar_seqlen+m+i*(n+m), L_Phi+m*m+i*(n+m)*(n+m), n+m, delta_z+m+i*(n+m), 1);
    }
    fwd_subst(tmp_optvar_seqlen+m+i*(n+m), L_Phi+m*m+i*(n+m)*(n+m), n, delta_z+m+i*(n+m), 1);
    
    bwd_subst(delta_z+m+i*(n+m), L_Phi_T+m*m+i*(n+m)*(n+m), n, tmp_optvar_seqlen+m+i*(n+m), 1);
    for (i = T-2; i > 0; i--){
        bwd_subst(delta_z+m+i*(n+m), L_Phi_T+m*m+i*(n+m)*(n+m), n+m, tmp_optvar_seqlen+m+i*(n+m), 1);
    }
    bwd_subst(delta_z+m+i*(n+m), L_Phi_T+m*m+i*(n+m)*(n+m), n+m, tmp_optvar_seqlen+m+i*(n+m), 1);
    bwd_subst(delta_z, L_Phi_T, m, tmp_optvar_seqlen, 1);
}

void form_delta_v(real_t delta_v[],
                  real_t *tmp_dual_seqlen,
                  const real_t L_Y[], const real_t L_Y_T[],
                  const real_t beta[],
                  const uint32_t T, const uint32_t n)
{
    uint32_t i;
    real_t tmp[n];
    mpcinc_mtx_scale(delta_v, beta, -1., T*n, 1);
    
    for (i = 0; i < T-1; i++) {
        fwd_subst(tmp_dual_seqlen+i*n, L_Y+2*i*n*n, n, delta_v+i*n, 1);
        mpcinc_mtx_multiply_mtx_vec(tmp, L_Y+2*i*n*n+n*n, tmp_dual_seqlen+i*n, n, n);
        mpcinc_mtx_substract_direct(delta_v+i*n+n, tmp, n, 1);
    }
    fwd_subst(tmp_dual_seqlen+i*n, L_Y+2*i*n*n, n, delta_v+i*n, 1); /*i=T-1*/
    
    for (i = T-1; i > 0; i--) {
        bwd_subst(delta_v+i*n, L_Y_T+2*i*n*n, n, tmp_dual_seqlen+i*n, 1);
        mpcinc_mtx_multiply_mtx_vec(tmp, L_Y_T+2*i*n*n-n*n, delta_v+i*n, n, n);
        mpcinc_mtx_substract_direct(tmp_dual_seqlen+i*n-n, tmp, n, 1);
    }
    bwd_subst(delta_v+i*n, L_Y_T+2*i*n*n, n, tmp_dual_seqlen+i*n, 1); /*i=0*/
}

void form_beta(real_t beta[],
               const real_t L_Phi[],
               const real_t L_Phi_T[],
               const real_t rd[], const real_t rp[],
               const uint32_t T,
               const real_t C[], const uint32_t n, const uint32_t m
               /*const real_t A[], const uint32_t n,
               const real_t B[], const uint32_t m*/)
{
    /* TODO beta lässt sich sicher auch parallel zu Y formen */
    real_t help1[T*(n+m)];
    real_t help2[T*(n+m)];
    real_t help3[T*n];
    
    fwd_subst(help1, L_Phi, T*(n+m), rd, 1);
    bwd_subst(help2, L_Phi_T, T*(n+m), help1, 1);
    mpcinc_mtx_multiply_mtx_vec(help3, C, help2, T*n, T*(n+m));
    mpcinc_mtx_substract(beta, help3, rp, T*n, 1);
}

void form_Y(real_t L_Y[], real_t *L_Y_T, real_t L_Phi[], real_t *L_Phi_T,
            const real_t Phi[],
            const uint32_t T,
            const real_t *A_B, const real_t *A_T_B_T, const uint32_t n,
            const real_t B[], const real_t *B_T, const uint32_t m,
            const real_t *eye_nm, const real_t *eye_n,
            real_t *PhiBlock, real_t *PhiBlock_I, real_t *last_PhiBlock_I,
            real_t *Qi_tilde, real_t *Y_bl, real_t *hilf1)
{
    uint32_t i, j, bl;
    /* hilf1 auch mehrmal als temporäre Variable für [n*n] und andere Größen verwendet */

    for (i = 0; i < T; i++){
        for (j = 0; j < (n+m)*(n+m); j++){  /* all blocks */
            last_PhiBlock_I[j] = PhiBlock_I[j];
        }
        if ( i == T-1){  /* last block i = T-1 */
            getBlock(PhiBlock, Phi, T*(n+m), m+i*(n+m), m+i*(n+m), n, n);
            bl = m*m + i*(n+m)*(n+m);
            cholesky(L_Phi+bl, PhiBlock, n);
            mpcinc_mtx_transpose(L_Phi_T+bl, L_Phi+bl, n, n);
            fwd_subst(hilf1, L_Phi+bl, n, eye_n, n);
            bwd_subst(Qi_tilde, L_Phi_T+bl, n, hilf1, n);
            form_Yii(Y_bl, A_B, n, n+m, last_PhiBlock_I, n+m, n+m, A_T_B_T, n+m, n, Qi_tilde, hilf1);
//             setBlock(Y, T*n, Y_bl, n, n, i*n, i*n);
            
            mpcinc_mtx_multiply_mtx_mtx(hilf1, L_Y+2*i*(n*n)-(n*n), L_Y_T+2*i*(n*n)-(n*n), n, n, n);
            mpcinc_mtx_scale_direct(hilf1, -1, n, n);
            mpcinc_mtx_add_direct(hilf1, Y_bl, n, n);
            cholesky(L_Y+2*i*(n*n), hilf1, n);
            mpcinc_mtx_transpose(L_Y_T+2*i*(n*n), L_Y+2*i*(n*n), n, n);
        }
        
        if (i < T-1){
            getBlock(PhiBlock, Phi, T*(n+m), m+i*(n+m), m+i*(n+m), n+m, n+m);
            bl = m*m + i*(n+m)*(n+m);
            cholesky(L_Phi+bl, PhiBlock, n+m);
            mpcinc_mtx_transpose(L_Phi_T+bl, L_Phi+bl, n+m, n+m);
            fwd_subst(hilf1, L_Phi+bl, n+m, eye_nm, n+m);
            bwd_subst(PhiBlock_I, L_Phi_T+bl, n+m, hilf1, n+m);
            getBlock(Qi_tilde, PhiBlock_I, n+m, 0, 0, n, n);
            
            if (i == 0){  /* first Block: i = 0 */
                getBlock(PhiBlock, Phi, T*(n+m), 0, 0, m, m);
                cholesky(L_Phi, PhiBlock, m);
                mpcinc_mtx_transpose(L_Phi_T, L_Phi, m, m);
                form_Y11(Y_bl, B, B_T, n, m, L_Phi, L_Phi_T, Qi_tilde, hilf1, hilf1+(m*n));
//                 setBlock(Y, T*n, Y_bl, n, n, i*n, i*n);
                /* hilf1 has size (n+m)*(n+m), so there is enough space for all */
                
                cholesky(L_Y+2*i*(n*n), Y_bl, n);
                mpcinc_mtx_transpose(L_Y_T+2*i*(n*n), L_Y+2*i*(n*n), n, n);
            }
            
            if (i > 0) {  /* not first block i != 0 */
                form_Yii(Y_bl, A_B, n, n+m, last_PhiBlock_I, n+m, n+m, A_T_B_T, n+m, n, Qi_tilde, hilf1);
//                 setBlock(Y, T*n, Y_bl, n, n, i*n, i*n);
                
                mpcinc_mtx_multiply_mtx_mtx(hilf1, L_Y+2*i*(n*n)-(n*n), L_Y_T+2*i*(n*n)-(n*n), n, n, n);
                mpcinc_mtx_scale_direct(hilf1, -1, n, n);
                mpcinc_mtx_add_direct(hilf1, Y_bl, n, n);
                cholesky(L_Y+2*i*(n*n), hilf1, n);
                mpcinc_mtx_transpose(L_Y_T+2*i*(n*n), L_Y+2*i*(n*n), n, n);
            }
            
            form_Y_i_ip1(Y_bl, A_T_B_T, n+m, n, PhiBlock_I);
//             setBlock(Y, T*n, Y_bl, n, n, i*n, (i+1)*n);
//             mpcinc_mtx_transpose(hilf1, Y_bl, n, n);
//             setBlock(Y, T*n, hilf1, n, n, (i+1)*n, i*n);
            
            fwd_subst(L_Y_T+2*i*(n*n)+(n*n), L_Y+2*i*(n*n), n, Y_bl, n);
            mpcinc_mtx_transpose(L_Y+2*i*(n*n)+(n*n), L_Y_T+2*i*(n*n)+(n*n), n, n);
        }
    }
}

void form_Yii(real_t sol[],
              const real_t A_B[], const uint32_t rowsA, const uint32_t colsA,
              const real_t last_PhiBlock_I[], const uint32_t rowsB, const uint32_t colsB,
              const real_t A_T_B_T[], const uint32_t rowsC, const uint32_t colsC,
              const real_t Qi[],
              real_t *tmp)
{
    mpcinc_mtx_multiply_mtx_mtx(tmp, last_PhiBlock_I, A_T_B_T, rowsB, colsB, colsC);
    mpcinc_mtx_multiply_mtx_mtx(sol, A_B, tmp, rowsA, colsA, colsC);
    mpcinc_mtx_add_direct(sol, Qi, rowsA, rowsA);
}

void form_Y_i_ip1(real_t sol[],
                  const real_t A_T_B_T[], const uint32_t rows,
                  const uint32_t cols,
                  const real_t QiSi_Block_I[])
{    
    mpcinc_mtx_multiply_mtx_mtx(sol, QiSi_Block_I, A_T_B_T, cols, rows, cols);
    mpcinc_mtx_scale_direct(sol, -1, cols, cols);
}

void form_Y11(real_t sol[],
              const real_t B[], const real_t *B_T,
              const uint32_t n, const uint32_t m,
              const real_t L_R0[], const real_t *L_R0_T,
              const real_t Q_tilde[],
              real_t *tmp1_mxn, real_t *tmp2_mxn)
{
    fwd_subst(tmp1_mxn, L_R0, m, B_T, n);
    bwd_subst(tmp2_mxn, L_R0_T, m, tmp1_mxn, n);
    mpcinc_mtx_multiply_mtx_mtx(sol, B, tmp2_mxn, n, m, n);
    mpcinc_mtx_add_direct(sol, Q_tilde, n, n);
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