/* Second order cone program solver based on a primal barrier interior point method */

#include <string.h>
#include "include/hhmpcipm.h"
/* static functions declaration */

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
        /* Solve system of linear equations to obtain the step direction */
        /* Find best step size (0...1] */
        /* Update z and x_k*/
    }
}

void hhmpc_ipm_check_valid(const struct hhmpc_ipm *ipm)
{
    /* TODO Über return FEHLER nachdenken, falls check_valid fehlschlägt.*/
}