/* Second order cone program solver based on a primal barrier interior point method */

#include <string.h>
#include "include/hhmpcipm.h"
#include "include/mpcincmtxops.h"
#include <hhmpcusefull.h>
/* static functions declaration */

static void residual(const struct hhmpc_ipm *ipm);
static void form_d(real_t *d, const real_t *P, const real_t *h, const real_t *z,
                   const uint32_t rowsP, const uint32_t colsP);

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
        /* Calculate the residual */
        residual(ipm);
        print_mtx(ipm->r_d, ipm->optvar_seqlen, 1);
        print_mtx(ipm->r_p, ipm->dual_seqlen, 1);
        /* Solve system of linear equations to obtain the step direction */
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
    real_t one[] = {1.};
    real_t *help = ipm->tmp1_optvar_seqlen;
    real_t *help2 = ipm->tmp2_dual_seqlen;
    
    form_d(ipm->d, ipm->P, ipm->h, ipm->z_opt, 36, ipm->optvar_seqlen);
    
    mpcinc_mtx_multiply_mtx_vec(help, ipm->P_T, ipm->d, ipm->optvar_seqlen, 36);
    mpcinc_mtx_scale(ipm->r_d, help, 0.015, ipm->optvar_seqlen, 1);
    mpcinc_mtx_mul_add(ipm->r_d, help, ipm->C_T, ipm->v_opt,
                       ipm->optvar_seqlen, ipm->dual_seqlen);
    mpcinc_mtx_mul_add(ipm->r_d, help, ipm->g, one,
                       ipm->optvar_seqlen, 1);
    
    mpcinc_mtx_scale(ipm->r_p, ipm->b, -1, ipm->dual_seqlen, 1);
    mpcinc_mtx_mul_add(ipm->r_p, help2, ipm->C, ipm->z_opt,
                       ipm->dual_seqlen, ipm->optvar_seqlen);
    
    print_mtx(ipm->P, 36, 15);
    print_mtx(ipm->d, 36, 1);
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
