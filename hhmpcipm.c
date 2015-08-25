/* Second order cone program solver based on a primal barrier interior point method */

#include <string.h>
#include "include/hhmpcipm.h"
#include "include/mpcincmtxops.h"
#include <hhmpcusefull.h>
/* static functions declaration */

static void residual(const struct hhmpc_ipm *ipm);

/* external functions definition */

void hhmpc_ipm_solve_problem(const struct hhmpc_ipm *ipm)
{
    uint32_t j;
    
    /*Check if initial value is valid*/
    hhmpc_ipm_check_valid(ipm);
    /*Take initial value*/
    memcpy(ipm->z_opt, ipm->z_ini, ipm->sizeof_optvar_seqlen);
    
    /*Improve z for a fixed number of steps j_in*/
    for (j = 0; j < *(ipm->j_in); j++) {
        /* Calculate the residual */
        residual(ipm);
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
    /*print_mtx(ipm->z_ini, ipm->optvar_seqlen, 1);*/
    real_t *help = ipm->tmp2_dual_seqlen;
    mpcinc_mtx_scale(ipm->r_p, ipm->b, -1, ipm->dual_seqlen, 1);
    mpcinc_mtx_mul_add(ipm->r_p, help, ipm->C, ipm->z_opt,
                       ipm->dual_seqlen, ipm->optvar_seqlen);
}
