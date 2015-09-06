#include "include/hhmpcsolve.h"


void solve_sysofleq(real_t delta_z[], real_t delta_v[],
                    const struct hhmpc_ipm *ipm,
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
    real_t *tmp1_optvar_seqlen = tmp_optvar_seqlen;
    real_t *tmp2_optvar_seqlen = ipm->tmp2_optvar_seqlen;
    real_t *tmp3_state_veclen = ipm->tmp3_state_veclen;
    
    real_t *PhiBlock = ipm->tmp_phibl1;
    real_t *PhiBlock_I = ipm->tmp_phibl2;
    real_t *PhiBlock_I_last = ipm->tmp_phibl3;
    real_t *Block_nxn1 = ipm->tmpQbl;
    real_t *Block_nxn2 = ipm->tmpYbl;
    real_t *tmp_optvar_veclenxoptvar_veclen = ipm->tmp10;
    real_t *L_Phi_blocks = ipm->tmp8_L_Phi;
    real_t *L_Phi_T_blocks = ipm->tmp9_L_Phi_T;
    
    form_Y(L_Y, L_Y_T, L_Phi_blocks, L_Phi_T_blocks,
           ipm->Phi, ipm->reg, T, ipm->A_B, ipm->A_B_T, n, B, B_T, m, eye_nm, eye_n,
           PhiBlock, PhiBlock_I, PhiBlock_I_last,
           Block_nxn1, Block_nxn2, tmp_optvar_veclenxoptvar_veclen);

    /* temporären Zeiger für beta sparen */
    form_beta(delta_v, L_Phi_blocks, L_Phi_T_blocks, rd, rp, T, C, n, m,
              tmp1_optvar_seqlen, tmp2_optvar_seqlen);
    form_delta_v(delta_v, tmp_dual_seqlen, tmp3_state_veclen, 
                 L_Y, L_Y_T, T, n);
    form_delta_z(delta_z, tmp1_optvar_seqlen, delta_v,
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
//     print_mtx(delta_z, T*(n+m), 1);
    
//     fwd_subst(tmp_optvar_seqlen, L_Phi, T*(n+m), delta_z, 1);
//     print_mtx(tmp_optvar_seqlen, T*(n+m), 1);
//     bwd_subst(delta_z, L_Phi_T, T*(n+m), tmp_optvar_seqlen, 1);
    
    fwd_subst(tmp_optvar_seqlen, L_Phi, m, delta_z, 1);
    for (i = 0; i < T-1; i++){
        fwd_subst(tmp_optvar_seqlen+m+i*(n+m), L_Phi+m*m+i*(n+m)*(n+m), n+m, delta_z+m+i*(n+m), 1);
    }
    fwd_subst(tmp_optvar_seqlen+m+i*(n+m), L_Phi+m*m+i*(n+m)*(n+m), n, delta_z+m+i*(n+m), 1);
//     print_mtx(tmp_optvar_seqlen, T*(n+m), 1);
    
    bwd_subst(delta_z+m+i*(n+m), L_Phi_T+m*m+i*(n+m)*(n+m), n, tmp_optvar_seqlen+m+i*(n+m), 1);
    for (i = T-2; i > 0; i--){
        bwd_subst(delta_z+m+i*(n+m), L_Phi_T+m*m+i*(n+m)*(n+m), n+m, tmp_optvar_seqlen+m+i*(n+m), 1);
    }
    bwd_subst(delta_z+m+i*(n+m), L_Phi_T+m*m+i*(n+m)*(n+m), n+m, tmp_optvar_seqlen+m+i*(n+m), 1);
    bwd_subst(delta_z, L_Phi_T, m, tmp_optvar_seqlen, 1);
//     print_mtx(delta_z, T*(n+m), 1);
}

void form_delta_v(real_t delta_v[],
                  real_t *tmp_dual_seqlen, real_t *tmp_n,
                  const real_t L_Y[], const real_t L_Y_T[],
                  const uint32_t T, const uint32_t n)
{
    uint32_t i;
    /* mpcinc_mtx_scale(delta_v, beta, -1., T*n, 1); */
    
    for (i = 0; i < T-1; i++) {
        fwd_subst(tmp_dual_seqlen+i*n, L_Y+2*i*n*n, n, delta_v+i*n, 1);
        mpcinc_mtx_multiply_mtx_vec(tmp_n, L_Y+2*i*n*n+n*n, tmp_dual_seqlen+i*n, n, n);
        mpcinc_mtx_substract_direct(delta_v+i*n+n, tmp_n, n, 1);
    }
    fwd_subst(tmp_dual_seqlen+i*n, L_Y+2*i*n*n, n, delta_v+i*n, 1); /*i=T-1*/
    
    for (i = T-1; i > 0; i--) {
        bwd_subst(delta_v+i*n, L_Y_T+2*i*n*n, n, tmp_dual_seqlen+i*n, 1);
        mpcinc_mtx_multiply_mtx_vec(tmp_n, L_Y_T+2*i*n*n-n*n, delta_v+i*n, n, n);
        mpcinc_mtx_substract_direct(tmp_dual_seqlen+i*n-n, tmp_n, n, 1);
    }
    bwd_subst(delta_v+i*n, L_Y_T+2*i*n*n, n, tmp_dual_seqlen+i*n, 1); /*i=0*/
}

void form_beta(real_t beta[],
               const real_t L_Phi[],
               const real_t L_Phi_T[],
               const real_t rd[], const real_t rp[],
               const uint32_t T,
               const real_t C[], const uint32_t n, const uint32_t m,
               real_t *tmp1, real_t *tmp2
               /*const real_t A[], const uint32_t n,
               const real_t B[], const uint32_t m*/)
{
    real_t beta_ref[] = {0.300426, 
69089400000.000000, 
59235000000.000572, 
60877400000.000000, 
59235000000.000000, 
51023000000.000000, 
-0.745263, 
-64100000000.000000, 
-64099999999.999336, 
-64100000000.000000, 
-64100000000.000000, 
-64100000000.000000, 
0.029387, 
-0.114214, 
-0.315186, 
-1.027924, 
-0.317260, 
-1987000000.000000, 
400.012495, 
0.000000, 
0.000833, 
0.000000, 
0.000000, 
0.000000, 
0.003654, 
1000000000.000000, 
999999999.999537, 
1000000000.000000, 
1000000000.000000, 
1000000000.000000, 
0.230973, 
52577033400.020332, 
45077835000.057045, 
46327701400.182976, 
45077835000.056473, 
39182189000.000000, 
4283542800.119049, 
30230535360.048935, 
25707990000.086334, 
27786644160.277542, 
25707990000.085663, 
19517046000.000000, 
0.318592, 
68398505999.901657, 
58642649999.728905, 
60268625999.114952, 
58642649999.726830, 
48801962999.999992, 
563078610062.702148, 
8082662632002.653320, 
7488166750000.762695, 
7761409192002.466797, 
7488166750000.761719, 
6608603500000.000000, 
0.003654, 
1000000000.000000, 
999999999.999537, 
1000000000.000000, 
1000000000.000000, 
1000000000.000000, 
0.230973, 
52577033400.020332, 
45077835000.057045, 
46327701400.182976, 
45077835000.056473, 
39182189000.000000, 
4283542800.119049, 
30230535360.048935, 
25707990000.086334, 
27786644160.277542, 
25707990000.085663, 
19517046000.000000, 
0.318592, 
68398505999.901657, 
58642649999.728905, 
60268625999.114952, 
58642649999.726830, 
48801962999.999992, 
563078610062.702148, 
8082662632002.653320, 
7488166750000.762695, 
7761409192002.466797, 
7488166750000.761719, 
6608603500000.000000, 
0.003654, 
1000000000.000000, 
999999999.999537, 
1000000000.000000, 
1000000000.000000, 
1000000000.000000, 
0.230973, 
52577033400.020332, 
45077835000.057045, 
46327701400.182976, 
45077835000.056473, 
39182189000.000000, 
4283542800.119049, 
30230535360.048935, 
25707990000.086334, 
27786644160.277542, 
25707990000.085663, 
19517046000.000000, 
0.318592, 
68398505999.901657, 
58642649999.728905, 
60268625999.114952, 
58642649999.726830, 
48801962999.999992, 
563078610062.702148, 
8082662632002.653320, 
7488166750000.762695, 
7761409192002.466797, 
7488166750000.761719, 
6608603500000.000000, 
0.003654, 
1000000000.000000, 
999999999.999537, 
1000000000.000000, 
1000000000.000000, 
1000000000.000000, 
-0.055984, 
-15512366599.979666, 
-13157164999.942953, 
-13549698599.817026, 
-13157164999.943523, 
-10840810999.999996, 
4283542800.533110, 
95330535360.048935, 
90807990000.086334, 
92886644160.277542, 
90807990000.085663, 
84617046000.000000, 
0.321184, 
68398506000.073349, 
58642650000.205833, 
60268626000.660202, 
58642650000.203758, 
51788962999.999992, 
563078610076.233887, 
8083662632002.653320, 
7489166750000.762695, 
7762409192002.466797, 
7489166750000.761719, 
6609603500000.000000, 
0.003654, 
1000000000.000000, 
999999999.999537, 
1000000000.000000, 
1000000000.000000, 
1000000000.000000};

    uint32_t i;
    /* TODO beta lässt sich sicher auch parallel zu Y formen */
    
//     fwd_subst(tmp1, L_Phi, T*(n+m), rd, 1);
//     bwd_subst(tmp2, L_Phi_T, T*(n+m), tmp1, 1);
    
    fwd_subst(tmp1, L_Phi, m, rd, 1);
    for (i = 0; i < T-1; i++){
        fwd_subst(tmp1+m+i*(n+m), L_Phi+m*m+i*(n+m)*(n+m), n+m, rd+m+i*(n+m), 1);
    }
    fwd_subst(tmp1+m+i*(n+m), L_Phi+m*m+i*(n+m)*(n+m), n, rd+m+i*(n+m), 1);
    
    bwd_subst(tmp2+m+i*(n+m), L_Phi_T+m*m+i*(n+m)*(n+m), n, tmp1+m+i*(n+m), 1);
    for (i = T-2; i > 0; i--){
        bwd_subst(tmp2+m+i*(n+m), L_Phi_T+m*m+i*(n+m)*(n+m), n+m, tmp1+m+i*(n+m), 1);
    }
    bwd_subst(tmp2+m+i*(n+m), L_Phi_T+m*m+i*(n+m)*(n+m), n+m, tmp1+m+i*(n+m), 1);
    bwd_subst(tmp2, L_Phi_T, m, tmp1, 1);
    
    mpcinc_mtx_multiply_mtx_vec(tmp1, C, tmp2, T*n, T*(n+m));
    mpcinc_mtx_substract(beta, tmp1, rp, T*n, 1);
//     print_mtx(beta, T*n, 1);
    printf("%d\n", mtx_cmp(beta, beta_ref, 150, 1e-1));
    printf("%f\t%f\n", beta[149], beta_ref[149]);
    mpcinc_mtx_scale_direct(beta, -1., T*n, 1);  /* return -beta */
}

void form_Y(real_t L_Y[], real_t *L_Y_T, real_t L_Phi[], real_t *L_Phi_T,
            real_t Phi[], const real_t *reg,
            const uint32_t T,
            const real_t *A_B, const real_t *A_T_B_T, const uint32_t n,
            const real_t B[], const real_t *B_T, const uint32_t m,
            const real_t *eye_nm, const real_t *eye_n,
            real_t *PhiBlock, real_t *PhiBlock_I, real_t *last_PhiBlock_I,
            real_t *Qi_tilde, real_t *Y_bl, real_t *hilf1)
{
    uint32_t i, j, ri, bl;
    /* hilf1 auch mehrmal als temporäre Variable für [n*n] und andere Größen verwendet */
    
    /* regularization (+epsilon*I) */
    for (ri = 0; ri < T*(n+m); ri++){
        (Phi+ri*T*(n+m)+ri)[0] += reg[0];
    }
    
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
            form_Yii(Y_bl, A_B, n, n+m, last_PhiBlock_I, n+m, n+m, A_T_B_T,
                     n+m, n, Qi_tilde, hilf1);
            /* regularization (-epsilon*I) */
            for (ri = 0; ri < n; ri++){
                (Y_bl+ri*n+ri)[0] += reg[0];
            }
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
                form_Y11(Y_bl, B, B_T, n, m, L_Phi, L_Phi_T, Qi_tilde,
                    hilf1, hilf1+(m*n));
                /* regularization (-epsilon*I) */
                for (ri = 0; ri < n; ri++){
                    (Y_bl+ri*n+ri)[0] += reg[0];
                }
//                 setBlock(Y, T*n, Y_bl, n, n, i*n, i*n);
                /* hilf1 has size (n+m)*(n+m), so there is enough space for all */
                
                cholesky(L_Y+2*i*(n*n), Y_bl, n);
                mpcinc_mtx_transpose(L_Y_T+2*i*(n*n), L_Y+2*i*(n*n), n, n);
            }
            
            if (i > 0) {  /* not first block i != 0 */
                form_Yii(Y_bl, A_B, n, n+m, last_PhiBlock_I, n+m, n+m, A_T_B_T,
                         n+m, n, Qi_tilde, hilf1);
                /* regularization (-epsilon*I) */
                for (ri = 0; ri < n; ri++){
                    (Y_bl+ri*n+ri)[0] += reg[0];
                }
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