/* Second order cone program solver based on a primal barrier interior point method */

#include <string.h>
#include "include/hhmpcipm.h"
#include "include/mpcincmtxops.h"
#include <hhmpcusefull.h>
#include "hhmpcsolve.h"
/* static functions declaration */

static void residual(const struct hhmpc_ipm *ipm,
                     const real_t *z, const real_t *v, const real_t *d,
                     const real_t kappa);
static void residual_norm(real_t *f, const real_t *r_d, const real_t* r_p,
                          const uint32_t optvar_seqlen, const uint32_t dual_seqlen);
static void bt_line_search(real_t *good_step, const struct hhmpc_ipm *ipm);
static void form_d(real_t *d, const real_t *P, const real_t *h, const real_t *z,
                   const uint32_t rowsP, const uint32_t colsP);
static void form_diag_d_sq(real_t *diag_d_sq, const real_t *d, const uint32_t dim);
static void form_Phi(real_t *Phi, real_t *help, real_t *tmp_Phi,
                     const real_t *H, const real_t *P_T, const real_t *P,
                     const struct hhmpc_ipm_P_hat *P_hat,
                     const real_t *d, const real_t *diag_d_sq,
                     const real_t kappa,
                     const uint32_t optvar_seqlen, const uint32_t nb_of_ueq_constr);
static void calc_kappa(real_t *kappa, const struct hhmpc_ipm *ipm,
                       const real_t *z);

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
    
    /* Update h for new xk */
    memcpy(ipm->P_of_z->h_hat, ipm->P_of_z->h, sizeof(real_t) * ipm->P_of_z->nb_lin_constr);
    
    print_mtx(ipm->h, 36, 1);
    /*Improve z for a fixed number of steps j_in*/
    for (j = 0; j < *(ipm->j_in); j++) {
        update(ipm->P_of_z, ipm->optvar_seqlen,
               t_solve_optvar_seqlen, t_optvar_seqlen);
        form_d(ipm->d, ipm->P, ipm->h, ipm->z_opt, ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
        form_diag_d_sq(ipm->diag_d_sq, ipm->d, ipm->nb_of_ueq_constr);
        form_Phi(ipm->Phi, ipm->tmp3_mtx_optvar_nb_of_ueq, tmp_Phi, ipm->H, ipm->P_T,
                 ipm->P, ipm->P_of_z , ipm->d, ipm->diag_d_sq, ipm->kappa[0],
                 ipm->optvar_seqlen, ipm->nb_of_ueq_constr);
        /* Calculate the residual */
        residual(ipm, ipm->z_opt, ipm->v_opt, ipm->d, ipm->kappa[0]);
        residual_norm(&f, ipm->r_d, ipm->r_p, ipm->optvar_seqlen, ipm->dual_seqlen);
        printf("res_norm = %f\n", f);
        
        /* Solve system of linear equations to obtain the step direction */
        solve_sysofleq(ipm->delta_z, ipm->delta_v, ipm, ipm->Phi, ipm->r_d, ipm->r_p,
                       ipm->C, ipm->C_T, ipm->A, ipm->A_T, ipm->B, ipm->B_T,
                       ipm->state_veclen, 2, ipm->horizon,
                       eye_nm, eye_n,
                       t_solve_optvar_seqlen,
                       t_solve_dual_seqlen,
                       t_L_Y, t_L_Y_T);
        /* Find best step size (0...1] */
        bt_line_search(ipm->st_size, ipm);
        print_mtx(ipm->st_size, 1,1);
        
        /* Update z */
        mpcinc_mtx_scale_direct(ipm->delta_z, ipm->st_size[0],
                                ipm->optvar_seqlen, 1);
        mpcinc_mtx_add_direct(ipm->z_opt, ipm->delta_z,
                              ipm->optvar_seqlen, 1);
        mpcinc_mtx_scale_direct(ipm->delta_v, ipm->st_size[0],
                                ipm->optvar_seqlen, 1);
        mpcinc_mtx_add_direct(ipm->v_opt, ipm->delta_v,
                              ipm->optvar_seqlen, 1);
        print_mtx(ipm->z_opt, ipm->optvar_seqlen, 1);
        print_mtx(ipm->v_opt, ipm->dual_seqlen, 1);
    }
    
    residual(ipm, ipm->z_opt, ipm->v_opt, ipm->d, ipm->kappa[0]);
    residual_norm(&f, ipm->r_d, ipm->r_p, ipm->optvar_seqlen, ipm->dual_seqlen);
    printf("res_norm = %.11f\n", f);
    /* Update x_k (und andere Parameter) */
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
}

void bt_line_search(real_t *st_size, const struct hhmpc_ipm *ipm)
{
    const real_t g_step = 1e-6;
    const real_t alpha = 0.4;
    const real_t beta = 0.6;
    real_t *help_z = ipm->tmp6_optvar_seqlen;
    real_t *help_v = ipm->tmp7_dual_seqlen;
    
    *st_size = 1.;
    
    real_t f_p;
    real_t f_p_g;
    real_t g_in_dir;
    
    
    residual_norm(&f_p, ipm->r_d, ipm->r_p,
                  ipm->optvar_seqlen, ipm->dual_seqlen);
    
    mpcinc_mtx_scale(help_z, ipm->delta_z, g_step, ipm->optvar_seqlen, 1);
    mpcinc_mtx_add_direct(help_z, ipm->z_opt, ipm->optvar_seqlen, 1);
    mpcinc_mtx_scale(help_v, ipm->delta_v, g_step, ipm->dual_seqlen, 1);
    mpcinc_mtx_add_direct(help_v, ipm->v_opt, ipm->dual_seqlen, 1);
    form_d(ipm->d, ipm->P, ipm->h, help_z, ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
    residual(ipm, help_z, help_v, ipm->d, ipm->kappa[0]);
    residual_norm(&f_p_g, ipm->r_d, ipm->r_p,
                  ipm->optvar_seqlen, ipm->dual_seqlen);
    g_in_dir = (f_p_g - f_p)/g_step;
    printf("%.8f\n", g_in_dir);
    
    mpcinc_mtx_scale(help_z, ipm->delta_z, st_size[0], ipm->optvar_seqlen, 1);
    mpcinc_mtx_add_direct(help_z, ipm->z_opt, ipm->optvar_seqlen, 1);
    hhmpc_ipm_check_valid(ipm, help_z);
    mpcinc_mtx_scale(help_v, ipm->delta_v, st_size[0], ipm->dual_seqlen, 1);
    mpcinc_mtx_add_direct(help_v, ipm->v_opt, ipm->dual_seqlen, 1);
    form_d(ipm->d, ipm->P, ipm->h, help_z, ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
    residual(ipm, help_z, help_v, ipm->d, ipm->kappa[0]);
    residual_norm(&f_p_g, ipm->r_d, ipm->r_p,
                  ipm->optvar_seqlen, ipm->dual_seqlen);
    while (hhmpc_ipm_check_valid(ipm, help_z) || (f_p_g > (f_p + alpha**st_size*g_in_dir)) )
    {
        *st_size *= beta;
        
        mpcinc_mtx_scale(help_z, ipm->delta_z, st_size[0], ipm->optvar_seqlen, 1);
        mpcinc_mtx_add_direct(help_z, ipm->z_opt, ipm->optvar_seqlen, 1);
        mpcinc_mtx_scale(help_v, ipm->delta_v, st_size[0], ipm->dual_seqlen, 1);
        mpcinc_mtx_add_direct(help_v, ipm->v_opt, ipm->dual_seqlen, 1);
        form_d(ipm->d, ipm->P, ipm->h, help_z, ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
        residual(ipm, help_z, help_v, ipm->d, ipm->kappa[0]);
        residual_norm(&f_p_g, ipm->r_d, ipm->r_p,
                      ipm->optvar_seqlen, ipm->dual_seqlen);
    }
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
    
    mpcinc_mtx_multiply_mtx_vec(help, ipm->P_T, d, ipm->optvar_seqlen, ipm->nb_of_ueq_constr);
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
    mpcinc_mtx_scale_direct(Phi, kappa, optvar, optvar);
    mpcinc_mtx_add_direct(Phi, H, optvar, optvar);  /* statt 2*H */
    mpcinc_mtx_add_direct(Phi, H, optvar, optvar);
    
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
    mpcinc_mtx_multiply_mtx_vec(kappa, ipm->g, z, 1, ipm->optvar_seqlen);
    kappa[0] *= 0.01/ipm->optvar_veclen;
}

