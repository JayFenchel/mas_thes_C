/* Second order cone program solver based on a primal barrier interior point method */

#include <string.h>
#include "include/hhmpcipm.h"
#include "include/mpcincmtxops.h"
#include <hhmpcusefull.h>
#include "hhmpcsolve.h"
/* static functions declaration */


static void residual_norm(real_t *f, const real_t *r_d, const real_t* r_p,
                          const uint32_t optvar_seqlen, const uint32_t dual_seqlen);
static void bt_line_search(real_t *good_step, const struct hhmpc_ipm *ipm);





/* external functions definition */

void hhmpc_ipm_solve_problem(const struct hhmpc_ipm *ipm)
{
    uint32_t j;
    real_t *t_solve_optvar_seqlen = ipm->tmp1_optvar_seqlen;
    real_t *t_optvar_seqlen = ipm->tmp2_optvar_seqlen;
    real_t *t_solve_dual_seqlen = ipm->tmp2_dual_seqlen;
    real_t *t_L_Y = ipm->tmp8_L_Y;
    real_t *t_L_Y_T = ipm->tmp9_L_Y_T;
    real_t *tmp_Phi = ipm->tmp4_mtx_optvar_optvar;
    real_t *eye_nm = ipm->eye_optvar_veclen;
    real_t *eye_n = ipm->eye_state_veclen;
    real_t f;
    real_t delta_ref1[] = {0.458037, 
-0.562717, 
0.000159, 
-0.037551, 
0.000023, 
0.000011, 
0.000011, 
-0.660641, 
-0.000011, 
-0.043997, 
-0.000023, 
-0.000023, 
-0.000011, 
-2.051853, 
0.000012, 
-0.136800, 
-0.000000, 
-0.000000, 
-0.000001, 
-0.824477, 
-0.000000, 
-0.054965, 
-0.000000, 
-0.000000, 
-0.000000, 
0.458032, 
0.000000, 
0.030536, 
0.000000, 
0.000000, 
0.000000, 
0.233336, 
-0.786311, 
0.000193, 
-0.052455, 
0.000011, 
0.000011, 
0.000011, 
-1.341861, 
0.034879, 
-0.089369, 
-0.000034, 
0.002296, 
-0.000034, 
-0.773407, 
-0.000059, 
-0.051530, 
0.000009, 
0.000005, 
-0.000003, 
-20.998529, 
4.578371, 
-1.392124, 
-0.001233, 
0.306444, 
0.000664, 
0.233334, 
0.000000, 
0.015555, 
0.000000, 
0.000000, 
0.000000, 
0.171631, 
-0.536494, 
0.000011, 
-0.035777, 
0.000000, 
0.000023, 
0.000017, 
-1.505526, 
0.083628, 
-0.100238, 
-0.000011, 
0.005548, 
-0.000068, 
-0.097955, 
-0.000011, 
-0.006487, 
0.000005, 
0.000006, 
-0.000007, 
-70.726167, 
13.222905, 
-4.696825, 
-0.001558, 
0.884995, 
-0.000581, 
0.171629, 
0.000000, 
0.011442, 
0.000000, 
0.000000, 
0.000000, 
0.087703, 
-0.253463, 
-0.000045, 
-0.016882, 
0.000011, 
0.000011, 
0.000006, 
-1.458717, 
0.116870, 
-0.097123, 
-0.000034, 
0.007765, 
-0.000080, 
0.124593, 
0.000031, 
0.008334, 
0.000011, 
0.000005, 
-0.000004, 
-141.388565, 
22.958106, 
-9.395985, 
-0.001458, 
1.534381, 
-0.002879, 
0.087702, 
0.000000, 
0.005847, 
0.000000, 
0.000000, 
0.000000, 
0.040631, 
-0.088374, 
-0.000002, 
-0.005881, 
-0.000001, 
0.000001, 
-0.000000, 
-1.389268, 
0.132591, 
-0.092468, 
-0.000005, 
0.008836, 
-0.000093, 
0.086212, 
0.000034, 
0.005756, 
0.000013, 
0.000007, 
-0.000001, 
-222.271886, 
32.518642, 
-14.776785, 
-0.002205, 
2.173635, 
-0.007107, 
0.040630, 
0.000000, 
0.002708, 
0.000000, 
0.000000, 
0.000000, 
66804.132534, 
11232.853757, 
-13.418792, 
-9.999824, 
-9.999589, 
-10.000000, 
-106860.302512, 
-10.000147, 
-9.999937, 
-10.000000, 
-10.000010, 
-10.000000, 
15103.581897, 
1216.206681, 
-1.000993, 
-9.999985, 
-9.999951, 
-10.000000, 
-1150.861084, 
-10.000001, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-442.153284, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
41201.242670, 
6441.616900, 
-10.230684, 
-9.999831, 
-9.999514, 
-10.000000, 
-54225.137249, 
-10.000100, 
-9.999957, 
-10.000000, 
-10.000007, 
-10.000000, 
5756.085217, 
549.863432, 
-6.602395, 
-10.000122, 
-10.000298, 
-10.000000, 
-855.870372, 
-10.000001, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-230.150731, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
20675.523931, 
3009.984360, 
-9.348347, 
-9.999685, 
-9.999525, 
-10.000000, 
-20614.267580, 
-10.000056, 
-9.999975, 
-10.000000, 
-10.000004, 
-10.000000, 
656.839909, 
150.402556, 
-9.609691, 
-10.000084, 
-10.000376, 
-10.000000, 
-575.788285, 
-10.000001, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-171.932011, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
6937.540500, 
891.132551, 
-9.643091, 
-9.999748, 
-9.999538, 
-10.000000, 
-2858.385358, 
-10.000021, 
-9.999991, 
-10.000000, 
-10.000001, 
-10.000000, 
-1129.854450, 
-10.006188, 
-10.572088, 
-10.000250, 
-10.000369, 
-10.000000, 
-332.454924, 
-10.000001, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-92.746833, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
169.346346, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
2720.693511, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-1003.472514, 
-10.005835, 
-10.360521, 
-10.000251, 
-10.000466, 
-10.000000, 
-141.341074, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-48.334457, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000};
    real_t delta_ref2[] = {0.458519, 
-0.563978, 
-0.000000, 
-0.037599, 
-0.000000, 
0.000000, 
0.000000, 
-0.660267, 
0.000000, 
-0.044018, 
-0.000000, 
-0.000000, 
0.000000, 
-2.054164, 
-0.000000, 
-0.136944, 
0.000000, 
0.000000, 
0.000000, 
-0.825334, 
0.000000, 
-0.055022, 
0.000000, 
0.000000, 
-0.000000, 
0.458519, 
0.000000, 
0.030568, 
-0.000000, 
0.000000, 
0.000000, 
0.233796, 
-0.788002, 
-0.000000, 
-0.052533, 
0.000000, 
-0.000000, 
0.000000, 
-1.341758, 
0.034967, 
-0.089451, 
0.000000, 
0.002331, 
-0.000000, 
-0.774599, 
0.000000, 
-0.051640, 
-0.000000, 
0.000000, 
0.000000, 
-20.920749, 
4.596421, 
-1.394717, 
-0.000000, 
0.306428, 
-0.000000, 
0.233796, 
0.000000, 
0.015586, 
0.000000, 
-0.000000, 
0.000000, 
0.172060, 
-0.537845, 
-0.000000, 
-0.035856, 
0.000000, 
0.000000, 
-0.000000, 
-1.505530, 
0.083823, 
-0.100369, 
0.000000, 
0.005588, 
0.000000, 
-0.098378, 
0.000000, 
-0.006559, 
-0.000000, 
0.000000, 
-0.000000, 
-70.562928, 
13.259996, 
-4.704195, 
0.000000, 
0.884000, 
-0.000000, 
0.172060, 
-0.000000, 
0.011471, 
0.000000, 
0.000000, 
-0.000000, 
0.087999, 
-0.254295, 
0.000000, 
-0.016953, 
-0.000000, 
-0.000000, 
0.000000, 
-1.458732, 
0.117169, 
-0.097249, 
0.000000, 
0.007811, 
-0.000000, 
0.124558, 
0.000000, 
0.008304, 
-0.000000, 
0.000000, 
0.000000, 
-141.161300, 
23.016473, 
-9.410753, 
0.000000, 
1.534432, 
-0.000000, 
0.087999, 
0.000000, 
0.005867, 
0.000000, 
0.000000, 
0.000000, 
0.040803, 
-0.088793, 
0.000000, 
-0.005920, 
0.000000, 
0.000000, 
0.000000, 
-1.389260, 
0.132935, 
-0.092617, 
0.000000, 
0.008862, 
-0.000000, 
0.086269, 
0.000000, 
0.005751, 
0.000000, 
0.000000, 
0.000000, 
-222.005530, 
32.599517, 
-14.800369, 
0.000000, 
2.173301, 
-0.000000, 
0.040803, 
0.000000, 
0.002720, 
-0.000000, 
0.000000, 
0.000000, 
66860.741470, 
11242.339122, 
-13.427048, 
-10.000000, 
-10.000000, 
-10.000000, 
-106962.949563, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
15119.833894, 
1217.542101, 
-0.991898, 
-10.000000, 
-10.000000, 
-10.000000, 
-1151.403354, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-442.612446, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
41245.802611, 
6448.510277, 
-10.235767, 
-10.000000, 
-10.000000, 
-10.000000, 
-54292.395863, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
5764.336271, 
550.627856, 
-6.595253, 
-10.000000, 
-10.000000, 
-10.000000, 
-856.413276, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-230.586993, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
20703.714129, 
3014.114985, 
-9.350336, 
-10.000000, 
-10.000000, 
-10.000000, 
-20650.249229, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
659.598689, 
150.686251, 
-9.605014, 
-10.000000, 
-10.000000, 
-10.000000, 
-576.273709, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-172.338956, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
6949.710220, 
892.731745, 
-9.643400, 
-10.000000, 
-10.000000, 
-10.000000, 
-2870.990908, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-1129.866349, 
-10.000000, 
-10.570138, 
-10.000000, 
-10.000000, 
-10.000000, 
-332.819713, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-93.026836, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
170.195568, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
2720.681465, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-1004.128454, 
-10.000000, 
-10.360202, 
-10.000000, 
-10.000000, 
-10.000000, 
-141.537914, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-48.497378, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000, 
-10.000000};
    
    /*Check if initial value is valid*/
    hhmpc_ipm_check_valid(ipm, ipm->z_ini);
    
    /*Take initial value*/
    memcpy(ipm->z_opt, ipm->z_ini, ipm->sizeof_optvar_seqlen);
    memcpy(ipm->v_opt, ipm->v_ini, ipm->sizeof_dual_seqlen);
    
    /* Calculate Kappa for every time_step */
    calc_kappa(ipm->kappa, ipm, ipm->z_opt);
    printf("calculated kappa = %.20f\n", ipm->kappa[0]);
    
    /* Update h for new xk */
    memcpy(ipm->P_of_z->h_hat, ipm->P_of_z->h, sizeof(real_t) * ipm->P_of_z->nb_lin_constr);
    /*
    print_mtx(ipm->h, ipm->nb_of_ueq_constr, 1);*/
    /*Improve z for a fixed number of steps j_in*/
    for (j = 0; j < *(ipm->j_in); j++) {
        update(ipm->P_of_z, ipm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
        form_d(ipm->d, ipm->P, ipm->h, ipm->z_opt,
               ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
        form_diag_d_sq(ipm->diag_d_sq, ipm->d, ipm->nb_of_ueq_constr);
        form_Phi(ipm->Phi, ipm->tmp3_mtx_optvar_nb_of_ueq, tmp_Phi, ipm->H,
                 ipm->P2_T, ipm->P2, ipm->P_of_z , ipm->d, ipm->diag_d_sq,
                 ipm->kappa[0], ipm->optvar_seqlen, ipm->nb_of_ueq_constr);
        /* Calculate the residual */
        residual(ipm, ipm->z_opt, ipm->v_opt, ipm->d, ipm->kappa[0]);
        residual_norm(&f, ipm->r_d, ipm->r_p, ipm->optvar_seqlen, ipm->dual_seqlen);
//         print_mtx(ipm->r_d, ipm->optvar_seqlen, 1);
        printf("res_norm = %f\n", f);
        
        /* Solve system of linear equations to obtain the step direction */
        solve_sysofleq(ipm->delta_z, ipm->delta_v, ipm, ipm->Phi, ipm->r_d, ipm->r_p,
                       ipm->C, ipm->C_T, ipm->A, ipm->A_T, ipm->B, ipm->B_T,
                       ipm->state_veclen, ipm->optvar_veclen-ipm->state_veclen, /* TODO m einführen */
                       ipm->horizon,
                       eye_nm, eye_n,
                       t_solve_optvar_seqlen,
                       t_solve_dual_seqlen,
                       t_L_Y, t_L_Y_T);
//         real_t tom1[ipm->optvar_seqlen*ipm->optvar_seqlen];
//         real_t tom2[ipm->optvar_seqlen];
//         cholesky(tom1, ipm->Phi, ipm->optvar_seqlen);
        
//         print_mtx(ipm->delta_v, ipm->dual_seqlen, 1);
        /* Find best step size (0...1] */
        bt_line_search(ipm->st_size, ipm);
        printf("st_size = %f\n", ipm->st_size[0]);
        
        /* Update z */
        mpcinc_mtx_scale_direct(ipm->delta_z, ipm->st_size[0],
                                ipm->optvar_seqlen, 1);
        mpcinc_mtx_add_direct(ipm->z_opt, ipm->delta_z,
                              ipm->optvar_seqlen, 1);
        mpcinc_mtx_scale_direct(ipm->delta_v, ipm->st_size[0],
                                ipm->optvar_seqlen, 1);
        mpcinc_mtx_add_direct(ipm->v_opt, ipm->delta_v,
                              ipm->optvar_seqlen, 1);
        /*
        print_mtx(ipm->z_opt, ipm->optvar_seqlen, 1);
        print_mtx(ipm->v_opt, ipm->dual_seqlen, 1);*/
    }
    ipm->kappa[0] = 95.;
    update(ipm->P_of_z, ipm->optvar_seqlen,
           t_solve_optvar_seqlen, t_optvar_seqlen);
    form_d(ipm->d, ipm->P, ipm->h, ipm->z_opt,
           ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
    residual(ipm, ipm->z_opt, ipm->v_opt, ipm->d, ipm->kappa[0]);
    residual_norm(&f, ipm->r_d, ipm->r_p, ipm->optvar_seqlen, ipm->dual_seqlen);
    printf("res_norm = %.11f\n", f);
    /* Update x_k (und andere Parameter) */
    memcpy(ipm->z_ini, ipm->z_opt, ipm->sizeof_optvar_seqlen);
    memcpy(ipm->v_ini, ipm->v_opt, ipm->sizeof_dual_seqlen);
}

uint32_t hhmpc_ipm_check_valid(const struct hhmpc_ipm *ipm, const real_t *z_check)
{
    uint32_t i;
    real_t *help1 = ipm->tmp4_nb_of_constr;
    real_t *help2 = ipm->tmp5_nb_of_constr;
    mpcinc_mtx_scale(help1, ipm->h, -1, ipm->nb_of_ueq_constr, 1);
    mpcinc_mtx_mul_add(help1, help2, ipm->P, z_check,
                       ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
//     printf("check:\n");
    for (i = 0; i < ipm->nb_of_ueq_constr; i++){
        /*printf("%f\n", help1[i]);*/
        if (help1[i] >= 0) {return 1;}
    }
    return 0;
}

void update(struct hhmpc_ipm_P_hat *P, const uint32_t optvar_seqlen,
            real_t *tmp1, real_t *tmp2)
{
    uint32_t i;
    struct hhmpc_ipm_qc *qc_i;
    struct hhmpc_ipm_socc *socc_i;
    /*  TODO P vorher 0 setzen */
    memcpy(P->P_hat, P->P, sizeof(real_t) * P->nb_lin_constr*optvar_seqlen);
    P->P_hat += P->nb_lin_constr*optvar_seqlen;  /* Pointer wird nicht nur lokal verändern */
    
    /* Determine rows for qc */
    for (i = 0; i < P->nb_qc; i++){
        qc_i = P->qc[i];
        mpcinc_mtx_multiply_mtx_mtx(P->P_hat+qc_i->par_0,
                                    qc_i->par, qc_i->Gamma,
                                    1, qc_i->dimGamma, qc_i->dimGamma);
        mpcinc_mtx_add_direct(P->P_hat+qc_i->par_0,
                              qc_i->beta, 1, qc_i->dimGamma);
        P->P_hat += optvar_seqlen;
    }
    /* Determine rows for socc */
    for (i = 0; i < P->nb_socc; i++){
        socc_i = P->socc[i];
        mpcinc_mtx_multiply_mtx_vec(tmp1, socc_i->A, socc_i->par,
                                    socc_i->rowsA, socc_i->colsA);
        mpcinc_mtx_add_direct(tmp1, socc_i->b, 1, socc_i->rowsA);
        mpcinc_mtx_add_direct(tmp1, socc_i->b, 1, socc_i->rowsA);
        mpcinc_mtx_multiply_mtx_mtx(P->P_hat+socc_i->par_0,
                                    tmp1, socc_i->A,
                                    1, socc_i->rowsA, socc_i->colsA);
        mpcinc_mtx_multiply_mtx_vec(tmp1, socc_i->c, socc_i->par,
                                    1, socc_i->colsA);
        tmp1[0] += 2*socc_i->d[0];
        mpcinc_mtx_scale(tmp2, socc_i->c, tmp1[0], socc_i->colsA, 1);
        mpcinc_mtx_substract_direct(P->P_hat+socc_i->par_0,
                                    tmp2, 1, socc_i->colsA);
        P->P_hat += optvar_seqlen;
    }

    P->P_hat -= (P->nb_lin_constr+P->nb_qc + P->nb_socc)*optvar_seqlen;
    mpcinc_mtx_transpose(P->P_hat_T, P->P_hat,
                         P->nb_lin_constr + P->nb_qc + P->nb_socc,
                         optvar_seqlen);
    
/**************************/
    memcpy(P->P2_hat, P->P, sizeof(real_t) * P->nb_lin_constr*optvar_seqlen);
    P->P2_hat += P->nb_lin_constr*optvar_seqlen;  /* Pointer wird nicht nur lokal verändern */
    
    /* Determine rows for qc */
    for (i = 0; i < P->nb_qc; i++){
        qc_i = P->qc[i];
        mpcinc_mtx_multiply_mtx_mtx(P->P2_hat+qc_i->par_0,
                                    qc_i->par, qc_i->Gamma,
                                    1, qc_i->dimGamma, qc_i->dimGamma);
        /* 2*z */
        mpcinc_mtx_scale_direct(P->P2_hat+qc_i->par_0, 2., 1, qc_i->dimGamma);
        mpcinc_mtx_add_direct(P->P2_hat+qc_i->par_0,
                              qc_i->beta, 1, qc_i->dimGamma);
        P->P2_hat += optvar_seqlen;
    }
    /* Determine rows for socc */
    for (i = 0; i < P->nb_socc; i++){
        socc_i = P->socc[i];
        mpcinc_mtx_multiply_mtx_vec(tmp1, socc_i->A, socc_i->par,
                                    socc_i->rowsA, socc_i->colsA);
        /* 2*z */
        mpcinc_mtx_scale_direct(tmp1, 2., 1, socc_i->rowsA);
        mpcinc_mtx_add_direct(tmp1, socc_i->b, 1, socc_i->rowsA);
        mpcinc_mtx_add_direct(tmp1, socc_i->b, 1, socc_i->rowsA);
        mpcinc_mtx_multiply_mtx_mtx(P->P2_hat+socc_i->par_0,
                                    tmp1, socc_i->A,
                                    1, socc_i->rowsA, socc_i->colsA);
        mpcinc_mtx_multiply_mtx_vec(tmp1, socc_i->c, socc_i->par,
                                    1, socc_i->colsA);
        /* 2*z */
        tmp1[0] *= 2.;
        tmp1[0] += 2*socc_i->d[0];
        mpcinc_mtx_scale(tmp2, socc_i->c, tmp1[0], socc_i->colsA, 1);
        mpcinc_mtx_substract_direct(P->P2_hat+socc_i->par_0,
                                    tmp2, 1, socc_i->colsA);
        P->P2_hat += optvar_seqlen;
    }
    
    P->P2_hat -= (P->nb_lin_constr+P->nb_qc + P->nb_socc)*optvar_seqlen;
    mpcinc_mtx_transpose(P->P2_hat_T, P->P2_hat,
                         P->nb_lin_constr + P->nb_qc + P->nb_socc,
                         optvar_seqlen);
}

void bt_line_search(real_t *st_size, const struct hhmpc_ipm *ipm)
{
    const real_t g_step = 1e-6;
    const real_t alpha = 0.4;
    const real_t beta = 0.6;
    real_t *help_z = ipm->tmp6_optvar_seqlen;
    real_t *help_v = ipm->tmp7_dual_seqlen;
    real_t *t_solve_optvar_seqlen = ipm->tmp1_optvar_seqlen;
    real_t *t_optvar_seqlen = ipm->tmp2_optvar_seqlen;
    
    *st_size = 1.;
    
    real_t f_p;
    real_t f_p_g;
    real_t g_in_dir;
    
    /*Save maintain value*/
    memcpy(help_z, ipm->z_opt, ipm->sizeof_optvar_seqlen);
    memcpy(help_v, ipm->v_opt, ipm->sizeof_dual_seqlen);
    
    
    residual_norm(&f_p, ipm->r_d, ipm->r_p,
                  ipm->optvar_seqlen, ipm->dual_seqlen);
    
    mpcinc_mtx_scale(ipm->z_opt, ipm->delta_z, g_step, ipm->optvar_seqlen, 1);
    mpcinc_mtx_add_direct(ipm->z_opt, help_z, ipm->optvar_seqlen, 1);
    mpcinc_mtx_scale(ipm->v_opt, ipm->delta_v, g_step, ipm->dual_seqlen, 1);
    mpcinc_mtx_add_direct(ipm->v_opt, help_v, ipm->dual_seqlen, 1);
    /* update P matrices */
    update(ipm->P_of_z, ipm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
    form_d(ipm->d, ipm->P, ipm->h, ipm->z_opt, ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
    residual(ipm, ipm->z_opt, ipm->v_opt, ipm->d, ipm->kappa[0]);
    residual_norm(&f_p_g, ipm->r_d, ipm->r_p,
                  ipm->optvar_seqlen, ipm->dual_seqlen);
    g_in_dir = (f_p_g - f_p)/g_step;
    printf("Grad in dir = %.8f\n", g_in_dir);
    
    mpcinc_mtx_scale(ipm->z_opt, ipm->delta_z, st_size[0], ipm->optvar_seqlen, 1);
    mpcinc_mtx_add_direct(ipm->z_opt, help_z, ipm->optvar_seqlen, 1);
    hhmpc_ipm_check_valid(ipm, ipm->z_opt);
    mpcinc_mtx_scale(ipm->v_opt, ipm->delta_v, st_size[0], ipm->dual_seqlen, 1);
    mpcinc_mtx_add_direct(ipm->v_opt, help_v, ipm->dual_seqlen, 1);
    /* update P matrices */
    update(ipm->P_of_z, ipm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
    form_d(ipm->d, ipm->P, ipm->h, ipm->z_opt, ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
    residual(ipm, ipm->z_opt, ipm->v_opt, ipm->d, ipm->kappa[0]);
    residual_norm(&f_p_g, ipm->r_d, ipm->r_p,
                  ipm->optvar_seqlen, ipm->dual_seqlen);
    while (hhmpc_ipm_check_valid(ipm, ipm->z_opt) || (f_p_g > (f_p + alpha**st_size*g_in_dir)) )
    {
        *st_size *= beta;
        
        mpcinc_mtx_scale(ipm->z_opt, ipm->delta_z, st_size[0], ipm->optvar_seqlen, 1);
        mpcinc_mtx_add_direct(ipm->z_opt, help_z, ipm->optvar_seqlen, 1);
        mpcinc_mtx_scale(ipm->v_opt, ipm->delta_v, st_size[0], ipm->dual_seqlen, 1);
        mpcinc_mtx_add_direct(ipm->v_opt, help_v, ipm->dual_seqlen, 1);
        /* update P matrices */
        update(ipm->P_of_z, ipm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
        form_d(ipm->d, ipm->P, ipm->h, ipm->z_opt, ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
        residual(ipm, ipm->z_opt, ipm->v_opt, ipm->d, ipm->kappa[0]);
        residual_norm(&f_p_g, ipm->r_d, ipm->r_p,
                      ipm->optvar_seqlen, ipm->dual_seqlen);
    }
    
    /*Load back maintain value*/
    memcpy(ipm->z_opt, help_z, ipm->sizeof_optvar_seqlen);
    memcpy(ipm->v_opt, help_v, ipm->sizeof_dual_seqlen);
}

void residual_norm(real_t *f, const real_t *r_d, const real_t *r_p,
                   const uint32_t ld, const uint32_t lp)
{
    uint32_t i;
    f[0] = 0.;
    for (i = 0; i < ld; i++){
        f[0] += r_d[i]*r_d[i];
    }
    for (i = 0; i < lp; i++){
        f[0] += r_p[i]*r_p[i];
    } 
}


void residual(const struct hhmpc_ipm *ipm,
              const real_t *z, const real_t *v, const real_t *d, const real_t kappa)
{
    real_t *help = ipm->tmp1_optvar_seqlen;
    real_t *help2 = ipm->tmp2_dual_seqlen;
    
    mpcinc_mtx_multiply_mtx_vec(help, ipm->P2_T, d, ipm->optvar_seqlen, ipm->nb_of_ueq_constr);
    mpcinc_mtx_scale(ipm->r_d, help, kappa, ipm->optvar_seqlen, 1);
    mpcinc_mtx_mul_add(ipm->r_d, help, ipm->C_T, v,
                       ipm->optvar_seqlen, ipm->dual_seqlen);
    mpcinc_mtx_add_direct(ipm->r_d, ipm->g,
                          ipm->optvar_seqlen, 1);
    mpcinc_mtx_multiply_mtx_vec(help, ipm->H, z,
                                ipm->optvar_seqlen, ipm->optvar_seqlen);
    mpcinc_mtx_scale_direct(help, 2, ipm->optvar_seqlen, 1);
    mpcinc_mtx_add_direct(ipm->r_d, help, ipm->optvar_seqlen, 1);
    
    mpcinc_mtx_scale(ipm->r_p, ipm->b, -1, ipm->dual_seqlen, 1);
    mpcinc_mtx_mul_add(ipm->r_p, help2, ipm->C, z,
                       ipm->dual_seqlen, ipm->optvar_seqlen);
}

void form_Phi(real_t *Phi, real_t *help, real_t *t_Phi,
              const real_t *H, const real_t *P_T, const real_t *P,
              const struct hhmpc_ipm_P_hat *P_hat,
              const real_t *d, const real_t *diag_d_sq,
              const real_t kappa,
              const uint32_t optvar, const uint32_t nb_of_ueq)
{    
    uint32_t i, j, pos_d, pos_Phi;
    struct hhmpc_ipm_qc *qc_i;
    struct hhmpc_ipm_socc *socc_i;

    mpcinc_mtx_multiply_mtx_mtx(help, P_T, diag_d_sq, optvar, nb_of_ueq, nb_of_ueq);
    mpcinc_mtx_multiply_mtx_mtx(Phi, help, P, optvar, nb_of_ueq, optvar);

    /* Additional terms for Phi resulting of the second derivative of qc */
    for (i = 0; i < P_hat->nb_qc; i++){
        qc_i = P_hat->qc[i];
        pos_d = P_hat->nb_lin_constr+i;
        mpcinc_mtx_scale(t_Phi, qc_i->Gamma, 2*d[pos_d],
                         qc_i->dimGamma, qc_i->dimGamma);
        
        for (j = 0; j < qc_i->dimGamma; j++){
            pos_Phi = (qc_i->par_0 + j)*optvar + qc_i->par_0;
            mpcinc_mtx_add_direct(Phi+pos_Phi, t_Phi+j*qc_i->dimGamma, 1, qc_i->dimGamma);
        }
    }
    /* Additional terms for Phi resulting of the second derivative of socc */
    for (i = 0; i < P_hat->nb_socc; i++){
        socc_i = P_hat->socc[i];
        pos_d = P_hat->nb_lin_constr + P_hat->nb_qc + i;
        mpcinc_mtx_scale(t_Phi, socc_i->AAmcc, 2*d[pos_d],
                         socc_i->colsA, socc_i->colsA);
        
        for (j = 0; j < socc_i->colsA; j++){
            pos_Phi = (socc_i->par_0 + j)*optvar + socc_i->par_0;
            mpcinc_mtx_add_direct(Phi+pos_Phi, t_Phi+j*socc_i->colsA, 1, socc_i->colsA);
        }
    }
    
    mpcinc_mtx_scale_direct(Phi, kappa, optvar, optvar);
    mpcinc_mtx_add_direct(Phi, H, optvar, optvar);  /* statt 2*H */
    mpcinc_mtx_add_direct(Phi, H, optvar, optvar);
}

void form_d(real_t *d, const real_t *P, const real_t *h, const real_t *z,
            const uint32_t rowsP, const uint32_t colsP)
{
    uint32_t i, j;
    for (i = 0; i < rowsP; i++){
        d[i] = h[i];
        for (j = 0; j < colsP; j++){
            d[i] -= P[i*colsP + j]*z[j];
        }
        d[i] = 1/d[i];
    }
}

void form_diag_d_sq(real_t *diag_d_sq, const real_t *d, const uint32_t dim)
{
    uint32_t i, j;
    for (i = 0; i < dim; i++){
        for (j = 0; j < dim; j++){
            diag_d_sq[i*dim+j] = (i == j) ? d[i]*d[i] : 0.;
        }
    }
}

void calc_kappa(real_t *kappa, const struct hhmpc_ipm *ipm, const real_t *z)
{
    real_t *tmp1 = ipm->tmp1_optvar_seqlen;
    real_t *tmp2 = ipm->tmp3_state_veclen; 
    
    mpcinc_mtx_multiply_mtx_mtx(tmp1, ipm->z_opt, ipm->H,
                                1, ipm->optvar_seqlen, ipm->optvar_seqlen);
    mpcinc_mtx_multiply_mtx_vec(tmp2, tmp1, ipm->z_opt, 1, ipm->optvar_seqlen);
    mpcinc_mtx_multiply_mtx_vec(kappa, ipm->g, z, 1, ipm->optvar_seqlen);
    kappa[0] += tmp2[0];
    kappa[0] *= 0.01/ipm->optvar_veclen;  /* TODO auf optvar_seqlen umstellen*/
}

