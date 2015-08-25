/* Second order cone program solver based on a primal barrier interior point method */

#include <string.h>
#include "include/hhmpcipm.h"
#include "include/mpcincmtxops.h"
#include <hhmpcusefull.h>
#include "hhmpcalg.h"
/* static functions declaration */

static void residual(const struct hhmpc_ipm *ipm);
static void form_d(real_t *d, const real_t *P, const real_t *h, const real_t *z,
                   const uint32_t rowsP, const uint32_t colsP);
static void form_diag_d_sq(real_t *diag_d_sq, const real_t *d, const uint32_t dim);
static void form_Phi(real_t *Phi, real_t *help,
                     const real_t *H, const real_t *P_T, const real_t *P, const real_t *diag_d_sq,
                     const uint32_t optvar_seqlen, const uint32_t nb_of_ueq_constr);

/* external functions definition */

void hhmpc_ipm_solve_problem(const struct hhmpc_ipm *ipm)
{
    uint32_t j;
    
    /*Check if initial value is valid*/
    hhmpc_ipm_check_valid(ipm);
    /*Take initial value*/
    memcpy(ipm->z_opt, ipm->z_ini, ipm->sizeof_optvar_seqlen);
    memcpy(ipm->v_opt, ipm->v_ini, ipm->sizeof_dual_seqlen);
    
    /*Improve z for a fixed number of steps j_in*/
    for (j = 0; j < *(ipm->j_in); j++) {
        form_d(ipm->d, ipm->P, ipm->h, ipm->z_opt, 36, ipm->optvar_seqlen);
        form_diag_d_sq(ipm->diag_d_sq, ipm->d, ipm->nb_of_ueq_constr);
        form_Phi(ipm->Phi, ipm->tmp3_mtx_optvar_nb_of_ueq, ipm->H, ipm->P_T,
                 ipm->P, ipm->diag_d_sq,
                 ipm->optvar_seqlen, ipm->nb_of_ueq_constr);
        print_mtx(ipm->Phi, ipm->optvar_seqlen, ipm->optvar_seqlen);
        /* Calculate the residual */
        residual(ipm);
        print_mtx(ipm->r_d, ipm->optvar_seqlen, 1);
        print_mtx(ipm->r_p, ipm->dual_seqlen, 1);
        /* Solve system of linear equations to obtain the step direction */
        
        solve_sysofleq(ipm->delta_z, ipm->delta_v, ipm->Phi, ipm->r_d, ipm->r_p,
                       ipm->C, ipm->A, ipm->B, ipm->state_veclen, 2, ipm->horizon);
        print_mtx(ipm->delta_z, ipm->optvar_seqlen, 1);
        print_mtx(ipm->delta_v, ipm->dual_seqlen, 1);
        /* Find best step size (0...1] */
        /* Update z */
    }
    /* Update x_k (und andere Parameter) */
}

void hhmpc_ipm_check_valid(const struct hhmpc_ipm *ipm)
{
    /* TODO Über return FEHLER nachdenken, falls check_valid fehlschlägt.*/
}

void residual(const struct hhmpc_ipm *ipm)
{
    real_t *help = ipm->tmp1_optvar_seqlen;
    real_t *help2 = ipm->tmp2_dual_seqlen;
    
    mpcinc_mtx_multiply_mtx_vec(help, ipm->P_T, ipm->d, ipm->optvar_seqlen, 36);
    mpcinc_mtx_scale(ipm->r_d, help, 0.015, ipm->optvar_seqlen, 1);
    mpcinc_mtx_mul_add(ipm->r_d, help, ipm->C_T, ipm->v_opt,
                       ipm->optvar_seqlen, ipm->dual_seqlen);
    mpcinc_mtx_add_direct(ipm->r_d, ipm->g,
                          ipm->optvar_seqlen, 1);
    mpcinc_mtx_multiply_mtx_vec(help, ipm->H, ipm->z_opt,
                                ipm->optvar_seqlen, ipm->optvar_seqlen);
    mpcinc_mtx_scale_direct(help, 2, ipm->optvar_seqlen, 1);
    mpcinc_mtx_add_direct(ipm->r_d, help, ipm->optvar_seqlen, 1);
    
    mpcinc_mtx_scale(ipm->r_p, ipm->b, -1, ipm->dual_seqlen, 1);
    mpcinc_mtx_mul_add(ipm->r_p, help2, ipm->C, ipm->z_opt,
                       ipm->dual_seqlen, ipm->optvar_seqlen);
    
    print_mtx(ipm->P, 36, 15);
    print_mtx(ipm->d, 36, 1);
}

void form_Phi(real_t *Phi, real_t *help,
              const real_t *H, const real_t *P_T, const real_t *P, const real_t *diag_d_sq,
              const uint32_t optvar, const uint32_t nb_of_ueq)
{
    mpcinc_mtx_multiply_mtx_mtx(help, P_T, diag_d_sq, optvar, nb_of_ueq, nb_of_ueq);
    mpcinc_mtx_multiply_mtx_mtx(Phi, help, P, optvar, nb_of_ueq, optvar);
    mpcinc_mtx_scale_direct(Phi, 0.015, optvar, optvar);
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
