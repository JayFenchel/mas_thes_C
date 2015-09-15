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

static void hhmpc_ipm_warm_start(const struct hhmpc_ipm *ipm);



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
        form_dsoft(ipm->dsoft, ipm->diag_d_soft,
                   ipm->roh, ipm->z_opt, ipm->hsoft,
                   ipm->Fusoft, ipm->Fxsoft, ipm->Ffsoft,
                   ipm->rowsFusoft, ipm->control_veclen, ipm->rowsFfsoft,
                   ipm->state_veclen, ipm->horizon);
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
//         print_mtx(ipm->delta_z, ipm->optvar_seqlen, 1);
        iterative_refinement(ipm);
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
                                ipm->dual_seqlen, 1);
        mpcinc_mtx_add_direct(ipm->v_opt, ipm->delta_v,
                              ipm->dual_seqlen, 1);
        /*
        print_mtx(ipm->z_opt, ipm->optvar_seqlen, 1);
        print_mtx(ipm->v_opt, ipm->dual_seqlen, 1);*/
    }
//     ipm->kappa[0] = 95.;
    update(ipm->P_of_z, ipm->optvar_seqlen,
           t_solve_optvar_seqlen, t_optvar_seqlen);
    form_d(ipm->d, ipm->P, ipm->h, ipm->z_opt,
           ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
    residual(ipm, ipm->z_opt, ipm->v_opt, ipm->d, ipm->kappa[0]);
    residual_norm(&f, ipm->r_d, ipm->r_p, ipm->optvar_seqlen, ipm->dual_seqlen);
    printf("res_norm = %.11f\n", f);
    /* Update x_k (und andere Parameter) */
    
    if (ipm->conf->warm_start) {
        hhmpc_ipm_warm_start(ipm);
    }
    /*
    memcpy(ipm->z_ini, ipm->z_opt, ipm->sizeof_optvar_seqlen);
    memcpy(ipm->v_ini, ipm->v_opt, ipm->sizeof_dual_seqlen);*/
}

void hhmpc_ipm_warm_start(const struct hhmpc_ipm *ipm)
{
    real_t *tmp = ipm->tmp3_state_veclen;
    mpcinc_mtx_shift_sequence(ipm->z_ini, ipm->z_opt, ipm->optvar_veclen,
            ipm->optvar_seqlen);
    mpcinc_mtx_multiply_mtx_vec((ipm->z_ini)+ipm->optvar_seqlen-ipm->state_veclen,
                                ipm->A, (ipm->z_opt)+ipm->optvar_seqlen-ipm->state_veclen,
                                ipm->state_veclen, ipm->state_veclen);
    mpcinc_mtx_mul_add((ipm->z_ini)+ipm->optvar_seqlen-ipm->state_veclen, tmp,
                                ipm->B, (ipm->z_opt)+ipm->optvar_seqlen-ipm->optvar_veclen,
                                ipm->state_veclen, ipm->control_veclen);    
    mpcinc_mtx_shift_sequence(ipm->v_ini, ipm->v_opt, ipm->state_veclen,
            ipm->dual_seqlen);
    
    return;
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

void iterative_refinement(const struct hhmpc_ipm *ipm)
{
    real_t *delta_rd = ipm->tmp1_optvar_seqlen;
    real_t *delta_rp = ipm->tmp2_dual_seqlen;
    real_t *delta_delta_v = ipm->tmp7_dual_seqlen;
    real_t *delta_delta_z = ipm->tmp4_mtx_optvar_optvar+ipm->dual_seqlen;
    real_t *tmp1 = ipm->tmp2_optvar_seqlen;
    real_t *tmp2 = ipm->tmp6_optvar_seqlen;
    real_t *tmp3 = ipm->tmp3_state_veclen;
    real_t *tmp4 = ipm->tmp4_mtx_optvar_optvar;
    
    real_t *L_Phi_blocks = ipm->tmp8_L_Phi;
    real_t *L_Phi_T_blocks = ipm->tmp9_L_Phi_T;
    real_t *L_Y = ipm->tmp8_L_Y;
    real_t *L_Y_T = ipm->tmp9_L_Y_T;
    
    mpcinc_mtx_multiply_mtx_vec(tmp1, ipm->C_T, ipm->delta_v,
                                ipm->optvar_seqlen, ipm->dual_seqlen);
    mpcinc_mtx_multiply_mtx_vec(delta_rd, ipm->Phi, ipm->delta_z,
                                ipm->optvar_seqlen, ipm->optvar_seqlen);
    mpcinc_mtx_add_direct(delta_rd, tmp1, ipm->optvar_seqlen, 1);
    mpcinc_mtx_add_direct(delta_rd, ipm->r_d, ipm->optvar_seqlen, 1);
    mpcinc_mtx_scale_direct(delta_rd, -1., ipm->optvar_seqlen, 1);
    
    mpcinc_mtx_multiply_mtx_vec(delta_rp, ipm->C, ipm->delta_z,
                                ipm->dual_seqlen, ipm->optvar_seqlen);
    mpcinc_mtx_add_direct(delta_rp, ipm->r_p, ipm->optvar_seqlen, 1);
    mpcinc_mtx_scale_direct(delta_rp, -1., ipm->dual_seqlen, 1);
   
//     print_mtx(delta_rp, ipm->dual_seqlen, 1);
        
    form_beta(delta_delta_v, L_Phi_blocks, L_Phi_T_blocks, delta_rd, delta_rp,
              ipm->horizon, ipm->C,
              ipm->state_veclen, ipm->optvar_veclen-ipm->state_veclen,
              tmp1, tmp2);
    form_delta_v(delta_delta_v, tmp4, tmp3, 
                 L_Y, L_Y_T, ipm->horizon, ipm->state_veclen);
    form_delta_z(delta_delta_z, tmp1, delta_delta_v,
                 L_Phi_blocks, L_Phi_T_blocks, delta_rd, ipm->C_T, ipm->horizon,
                 ipm->state_veclen, ipm->optvar_veclen-ipm->state_veclen);
    
    
//     print_mtx(delta_delta_z, ipm->optvar_seqlen, 1);
    mpcinc_mtx_substract_direct(ipm->delta_z, delta_delta_z, ipm->optvar_seqlen, 1);
    mpcinc_mtx_substract_direct(ipm->delta_v, delta_delta_v, ipm->dual_seqlen, 1);
//     print_mtx(ipm->delta_z, ipm->optvar_seqlen, 1);
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
    
    mpcinc_mtx_scale(ipm->z_opt, ipm->delta_z, st_size[0],
                     ipm->optvar_seqlen, 1);
    mpcinc_mtx_add_direct(ipm->z_opt, help_z,
                          ipm->optvar_seqlen, 1);
    mpcinc_mtx_scale(ipm->v_opt, ipm->delta_v, st_size[0],
                     ipm->dual_seqlen, 1);
    mpcinc_mtx_add_direct(ipm->v_opt, help_v,
                          ipm->dual_seqlen, 1);
    hhmpc_ipm_check_valid(ipm, ipm->z_opt);
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
        
        mpcinc_mtx_scale(ipm->z_opt, ipm->delta_z, st_size[0],
                         ipm->optvar_seqlen, 1);
        mpcinc_mtx_add_direct(ipm->z_opt, help_z,
                              ipm->optvar_seqlen, 1);
        mpcinc_mtx_scale(ipm->v_opt, ipm->delta_v, st_size[0],
                         ipm->dual_seqlen, 1);
        mpcinc_mtx_add_direct(ipm->v_opt, help_v,
                              ipm->dual_seqlen, 1);
        /* update P matrices */
        update(ipm->P_of_z, ipm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
        form_d(ipm->d, ipm->P, ipm->h, ipm->z_opt,
               ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
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
    real_t *tmp1_os = ipm->tmp1_res_os;
    real_t *tmp2_os = ipm->tmp2_res_os;
    real_t *tmp3_ds = ipm->tmp3_res_ds;
    
    mpcinc_mtx_multiply_mtx_vec(tmp1_os, ipm->P2_T, d, ipm->optvar_seqlen, ipm->nb_of_ueq_constr);
    mpcinc_mtx_scale(ipm->r_d, tmp1_os, kappa, ipm->optvar_seqlen, 1);
    mpcinc_mtx_mul_add(ipm->r_d, tmp1_os, ipm->C_T, v,
                       ipm->optvar_seqlen, ipm->dual_seqlen);
    mpcinc_mtx_add_direct(ipm->r_d, ipm->g,
                          ipm->optvar_seqlen, 1);
    mpcinc_mtx_substract(tmp2_os, z, ipm->zref, ipm->optvar_seqlen, 1);
    mpcinc_mtx_multiply_mtx_vec(tmp1_os, ipm->H, tmp2_os,
                                ipm->optvar_seqlen, ipm->optvar_seqlen);
    mpcinc_mtx_scale_direct(tmp1_os, 2, ipm->optvar_seqlen, 1);
    mpcinc_mtx_add_direct(ipm->r_d, tmp1_os, ipm->optvar_seqlen, 1);
    
    mpcinc_mtx_scale(ipm->r_p, ipm->b, -1, ipm->dual_seqlen, 1);
    mpcinc_mtx_mul_add(ipm->r_p, tmp3_ds, ipm->C, z,
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

void form_dsoft(real_t *ds, real_t *diags,
                const real_t *roh, const real_t *z, const real_t *hs,
                const real_t *Fus, const real_t *Fxs, const real_t *Ffs,
                const uint32_t rowsFus, const uint32_t c_veclen,
                const uint32_t rowsFfs, const uint32_t s_veclen,
                const uint32_t T)
{
    real_t rd_soft[T*(c_veclen+s_veclen)];
    uint32_t k, i, j;
    for (i = 0; i < rowsFus; i++){
        ds[i] = hs[i];
        for (j = 0; j < c_veclen; j++){
            ds[i] -= (Fus+i*c_veclen)[j]*z[j];
        }
        ds[i] = smpl_pow(E, roh[0]*ds[i]);
        diags[i] = ds[i];
        ds[i] = 1 / (1 + ds[i]);
        diags[i] = diags[i] * ds[i] * ds[i];
    }
    for (k = 1; k < T; k++){
        for (i = k*rowsFus; i < (k+1)*rowsFus; i++){
            ds[i] = hs[i];
            for (j = 0; j < s_veclen; j++){
                ds[i] -= (Fxs+(i-k*rowsFus)*s_veclen)[j]*(z+k*c_veclen+(k-1)*s_veclen)[j];
            }
            for (j = 0; j < c_veclen; j++){
                ds[i] -= (Fus+(i-k*rowsFus)*c_veclen)[j]*(z+k*c_veclen+(k)*s_veclen)[j];
            }
            ds[i] = smpl_pow(E, roh[0]*ds[i]);
            diags[i] = ds[i];
            ds[i] = 1 / (1 + ds[i]);
            diags[i] = diags[i] * ds[i] * ds[i];
        }
    }
    for (i = T*rowsFus; i < T*rowsFus + rowsFfs; i++){
        ds[i] = hs[i];
        for (j = 0; j < s_veclen; j++){
            ds[i] -= (Ffs+(i-T*rowsFus)*s_veclen)[j]*(z+T*c_veclen+(T-1)*s_veclen)[j];
        }
        ds[i] = smpl_pow(E, roh[0]*ds[i]);
        diags[i] = ds[i];
        ds[i] = 1 / (1 + ds[i]);
        diags[i] = diags[i] * ds[i] * ds[i];
    }
//     print_mtx(ds, T*rowsFus+rowsFfs, 1);
//     print_mtx(diags, T*rowsFus+rowsFfs, 1);
    
    for (i = 0; i < c_veclen; i++){
        rd_soft[i] = 0;
        for (j = 0; j < rowsFus; j++){
            rd_soft[i] += Fus[i+j*c_veclen]*ds[j];
        }
    }
    
    for (i = 0; i < c_veclen; i++){
        rd_soft[i] = 0;
        for (j = 0; j < rowsFus; j++){
            rd_soft[i] += Fus[i+j*c_veclen]*ds[j];
        }
    }
    for (i = T*c_veclen+(T-1)*s_veclen; i < T*(c_veclen+s_veclen); i++){
        rd_soft[i] = 0.;
        for (j = 0; j < rowsFfs; j++){
            rd_soft[i] += Ffs[i-(T*c_veclen+(T-1)*s_veclen)+j*s_veclen] * (ds+T*rowsFus)[j];
        }
    }
    print_mtx(rd_soft, T*(c_veclen + s_veclen), 1);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
    real_t *tmp3 = ipm->tmp2_optvar_seqlen;
    real_t *tmp2 = ipm->tmp3_state_veclen; 
    
    mpcinc_mtx_substract(tmp3, ipm->z_opt, ipm->zref, ipm->optvar_seqlen, 1);
    
    mpcinc_mtx_multiply_mtx_mtx(tmp1, ipm->z_opt, ipm->H,
                                1, ipm->optvar_seqlen, ipm->optvar_seqlen);
    mpcinc_mtx_multiply_mtx_vec(tmp2, tmp1, ipm->z_opt, 1, ipm->optvar_seqlen);
    mpcinc_mtx_multiply_mtx_vec(kappa, ipm->g, z, 1, ipm->optvar_seqlen);
    kappa[0] += tmp2[0];
    kappa[0] *= 0.01/ipm->optvar_veclen;  /* TODO auf optvar_seqlen umstellen*/
}

