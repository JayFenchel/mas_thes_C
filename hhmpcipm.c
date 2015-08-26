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
static void form_Phi(real_t *Phi, real_t *help,
                     const real_t *H, const real_t *P_T, const real_t *P, const real_t *diag_d_sq,
                     const real_t kappa,
                     const uint32_t optvar_seqlen, const uint32_t nb_of_ueq_constr);
static void calc_kappa(real_t *kappa, const struct hhmpc_ipm *ipm,
                       const real_t *z);

/* external functions definition */

void hhmpc_ipm_solve_problem(const struct hhmpc_ipm *ipm)
{
    uint32_t j;
    real_t f;
    
    /*Check if initial value is valid*/
    hhmpc_ipm_check_valid(ipm, ipm->z_ini);
    
    /*Take initial value*/
    memcpy(ipm->z_opt, ipm->z_ini, ipm->sizeof_optvar_seqlen);
    memcpy(ipm->v_opt, ipm->v_ini, ipm->sizeof_dual_seqlen);
    
    /* Calculate Kappa once for every time_step */
    calc_kappa(ipm->kappa, ipm, ipm->z_opt);
    
    /*Improve z for a fixed number of steps j_in*/
    for (j = 0; j < *(ipm->j_in); j++) {
        form_d(ipm->d, ipm->P, ipm->h, ipm->z_opt, ipm->nb_of_ueq_constr, ipm->optvar_seqlen);
        form_diag_d_sq(ipm->diag_d_sq, ipm->d, ipm->nb_of_ueq_constr);
        form_Phi(ipm->Phi, ipm->tmp3_mtx_optvar_nb_of_ueq, ipm->H, ipm->P_T,
                 ipm->P, ipm->diag_d_sq, ipm->kappa[0],
                 ipm->optvar_seqlen, ipm->nb_of_ueq_constr);
        /* Calculate the residual */
        residual(ipm, ipm->z_opt, ipm->v_opt, ipm->d, ipm->kappa[0]);
        residual_norm(&f, ipm->r_d, ipm->r_p, ipm->optvar_seqlen, ipm->dual_seqlen);
        printf("res_norm = %f\n", f);
        
        /* Solve system of linear equations to obtain the step direction */
        solve_sysofleq(ipm->delta_z, ipm->delta_v, ipm->Phi, ipm->r_d, ipm->r_p,
                       ipm->C, ipm->A, ipm->B, ipm->state_veclen, 2, ipm->horizon);
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
    printf("res_norm = %f\n", f);
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
    printf("check:\n");
    for (i = 0; i < ipm->nb_of_ueq_constr; i++){
        /*printf("%f\n", help1[i]);*/
        if (help1[i] >= 0) {return 1;}
    }
    return 0;
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

void form_Phi(real_t *Phi, real_t *help,
              const real_t *H, const real_t *P_T, const real_t *P, const real_t *diag_d_sq,
              const real_t kappa,
              const uint32_t optvar, const uint32_t nb_of_ueq)
{
    mpcinc_mtx_multiply_mtx_mtx(help, P_T, diag_d_sq, optvar, nb_of_ueq, nb_of_ueq);
    mpcinc_mtx_multiply_mtx_mtx(Phi, help, P, optvar, nb_of_ueq, optvar);
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
    mpcinc_mtx_multiply_mtx_vec(kappa, ipm->g, z, 1, ipm->optvar_seqlen);
    kappa[0] *= 0.01/ipm->optvar_veclen;
}

