#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */

#include "mc04types.h"
#include "include/cjson.h"
#include "include/hhmpcipmdynmem.h"
#include "include/mpcincmtxops.h"
#include "include/hhmpcusefull.h"

/* Static functions declarations */
hhmpc_dynmem_error_t hhmpc_ipm_parse_elements(struct hhmpc_ipm *ipm, cJSON *data);

/* Extern function definitions */

struct hhmpc_ipm *hhmpc_ipm_allocate_solver(void)
{
    struct hhmpc_ipm *ipm = (struct hhmpc_ipm*)malloc(sizeof(struct hhmpc_ipm));
    if (NULL == ipm) {return NULL;}
    struct hhmpc_ipm_conf *conf = (struct hhmpc_ipm_conf*)malloc(sizeof(struct hhmpc_ipm_conf));
    if (NULL == conf) {return NULL;}
    ipm->conf = conf;
    struct hhmpc_ipm_P_hat *P = (struct hhmpc_ipm_P_hat*)malloc(sizeof(struct hhmpc_ipm_P_hat));
    if (NULL == P) {return NULL;}
    ipm->P_of_z = P;
    
    ipm->kappa = (real_t *)malloc(sizeof(real_t));
    if (NULL == ipm->kappa) {return NULL;}
    
    return ipm;
}

hhmpc_dynmem_error_t hhmpc_ipm_setup_solver(struct hhmpc_ipm *ipm,
                                            struct hhmpc_socp_prb *prb,
                                            char *fname)
{
    hhmpc_dynmem_error_t ret;
    cJSON *data;
    uint32_t i;
    real_t *tmp;
    
    data = hhmpc_dynmem_get_data(fname);
    if (NULL == data) {return HHMPC_DYNMEM_FAIL;}
    ret = hhmpc_ipm_parse_elements(ipm, data);
    if (HHMPC_DYNMEM_OK != ret) {return ret;}
    
    ipm->b = prb->b->data;
    /*ipm->h = prb->h->data;*/
    ipm->g = prb->g->data;
    ipm->A = prb->A->data;
    ipm->B = prb->B->data;
    ipm->C = prb->C->data;
    ipm->H = prb->H->data;
    ipm->z_ini = prb->z_ini->data;
    ipm->v_ini = prb->v_ini->data;
    ipm->zref = prb->zref->data;

    ipm->optvar_seqlen = ipm->optvar_veclen * ipm->horizon;
    ipm->dual_seqlen = ipm->state_veclen * ipm->horizon;
    ipm->control_veclen = ipm->optvar_veclen - ipm->state_veclen;
    ipm->sizeof_optvar_seqlen = sizeof(real_t) * ipm->optvar_seqlen;
    ipm->sizeof_dual_seqlen = sizeof(real_t) * ipm->dual_seqlen;
/*    
    ipm->z_ini = (real_t *)malloc(ipm->sizeof_optvar_seqlen);
    if (NULL == ipm->z_ini) {return HHMPC_DYNMEM_FAIL;}*/
    ipm->z_opt = (real_t *)malloc(ipm->sizeof_optvar_seqlen);
    if (NULL == ipm->z_opt) {return HHMPC_DYNMEM_FAIL;}
    ipm->delta_z = (real_t *)malloc(ipm->sizeof_optvar_seqlen);
    if (NULL == ipm->delta_z) {return HHMPC_DYNMEM_FAIL;}
/*    
    ipm->v_ini = (real_t *)malloc(ipm->sizeof_dual_seqlen);
    if (NULL == ipm->v_ini) {return HHMPC_DYNMEM_FAIL;}*/
    ipm->v_opt = (real_t *)malloc(ipm->sizeof_dual_seqlen);
    if (NULL == ipm->v_opt) {return HHMPC_DYNMEM_FAIL;}
    ipm->delta_v = (real_t *)malloc(ipm->sizeof_dual_seqlen);
    if (NULL == ipm->delta_v) {return HHMPC_DYNMEM_FAIL;}
    
    /* u_k points on first */
    prb->u_k->data = ipm->z_opt;
    prb->u_k->rows = ipm->control_veclen;
    prb->u_k->cols = 1;
    
    ipm->P_of_z->h = prb->h->data;
    ipm->P_of_z->P = prb->P->data;
    ipm->P_of_z->nb_lin_constr = prb->P->rows;
    ipm->P_of_z->nb_socc = prb->nb_socc;
    ipm->P_of_z->nb_qc = prb->nb_qc;
    ipm->nb_of_ueq_constr = prb->P->rows+prb->nb_socc+prb->nb_qc;
        
    ipm->P_of_z->h_hat =
            (real_t*)malloc(sizeof(real_t) * ipm->nb_of_ueq_constr);
    if (NULL == ipm->P_of_z->h_hat) {return HHMPC_DYNMEM_FAIL;}
    ipm->P_of_z->socc = (struct hhmpc_ipm_socc**)calloc(prb->nb_socc, sizeof(struct hhmpc_ipm_socc*));
    if (NULL == ipm->P_of_z->socc) {return HHMPC_DYNMEM_FAIL;}
    for (i = 0; i < prb->nb_socc; i++){
        ipm->P_of_z->socc[i] = (struct hhmpc_ipm_socc*)malloc(sizeof(struct hhmpc_ipm_socc));
        if (NULL == ipm->P_of_z->socc[i]) {return HHMPC_DYNMEM_FAIL;}
        ipm->P_of_z->socc[i]->rowsA = prb->socc[i]->A->rows;
        ipm->P_of_z->socc[i]->colsA = prb->socc[i]->A->cols;
        ipm->P_of_z->socc[i]->A = prb->socc[i]->A->data;
        ipm->P_of_z->socc[i]->A_T =
                (real_t *)malloc(sizeof(real_t) * prb->socc[i]->A->cols*prb->socc[i]->A->rows);
        if (NULL == ipm->P_of_z->socc[i]->A_T) {return HHMPC_DYNMEM_FAIL;}
        mpcinc_mtx_transpose(ipm->P_of_z->socc[i]->A_T, ipm->P_of_z->socc[i]->A,
                             prb->socc[i]->A->rows, prb->socc[i]->A->cols);
        ipm->P_of_z->socc[i]->b = prb->socc[i]->b->data;
        ipm->P_of_z->socc[i]->c = prb->socc[i]->c->data;
        ipm->P_of_z->socc[i]->d = prb->socc[i]->d->data;
        ipm->P_of_z->socc[i]->AAmcc=
                (real_t *)malloc(sizeof(real_t) * prb->socc[i]->A->cols*prb->socc[i]->A->cols);
        if (NULL == ipm->P_of_z->socc[i]->AAmcc) {return HHMPC_DYNMEM_FAIL;}
        mpcinc_mtx_multiply_mtx_mtx(ipm->P_of_z->socc[i]->AAmcc,
                                    ipm->P_of_z->socc[i]->A_T,
                                    ipm->P_of_z->socc[i]->A,
                                    ipm->P_of_z->socc[i]->colsA,
                                    ipm->P_of_z->socc[i]->rowsA,
                                    ipm->P_of_z->socc[i]->colsA);
        real_t cc[ipm->P_of_z->socc[i]->colsA*ipm->P_of_z->socc[i]->colsA];
        mpcinc_mtx_multiply_mtx_mtx(cc,
                                    ipm->P_of_z->socc[i]->c,
                                    ipm->P_of_z->socc[i]->c,
                                    ipm->P_of_z->socc[i]->colsA, 1,
                                    ipm->P_of_z->socc[i]->colsA);
        mpcinc_mtx_substract_direct(ipm->P_of_z->socc[i]->AAmcc, cc,
                                    ipm->P_of_z->socc[i]->colsA,
                                    ipm->P_of_z->socc[i]->colsA);
        
        ipm->P_of_z->socc[i]->par = ipm->z_opt+prb->socc[i]->par_0;
        ipm->P_of_z->socc[i]->par_0 = prb->socc[i]->par_0;
        ipm->P_of_z->socc[i]->par_l = prb->socc[i]->par_l;
        tmp = ipm->P_of_z->h_hat+prb->P->rows+prb->nb_qc+i;
        mpcinc_mtx_multiply_mtx_vec(tmp, ipm->P_of_z->socc[i]->b, ipm->P_of_z->socc[i]->b,
                                    1, ipm->P_of_z->socc[i]->rowsA);
        tmp[0] = ipm->P_of_z->socc[i]->d[0]*ipm->P_of_z->socc[i]->d[0] - tmp[0];
    }
    ipm->P_of_z->qc = (struct hhmpc_ipm_qc**)calloc(prb->nb_qc, sizeof(struct hhmpc_ipm_qc*));
    if (NULL == ipm->P_of_z->qc) {return HHMPC_DYNMEM_FAIL;}
    for (i = 0; i < prb->nb_qc; i++){
        ipm->P_of_z->qc[i] = (struct hhmpc_ipm_qc*)malloc(sizeof(struct hhmpc_ipm_qc));
        if (NULL == ipm->P_of_z->qc[i]) {return HHMPC_DYNMEM_FAIL;}
        ipm->P_of_z->qc[i]->dimGamma = prb->qc[i]->Gamma->rows;
        ipm->P_of_z->qc[i]->Gamma = prb->qc[i]->Gamma->data;
        ipm->P_of_z->qc[i]->beta = prb->qc[i]->beta->data;
        ipm->P_of_z->qc[i]->alpha = prb->qc[i]->alpha->data;
        ipm->P_of_z->qc[i]->par = ipm->z_opt+prb->qc[i]->par_0;
        ipm->P_of_z->qc[i]->par_0 = prb->qc[i]->par_0;
        ipm->P_of_z->qc[i]->par_l = prb->qc[i]->par_l;
        tmp = ipm->P_of_z->h_hat+prb->P->rows+i;
        tmp[0] = ipm->P_of_z->qc[i]->alpha[0];
    }
    ipm->P_of_z->P_hat =
            (real_t *)malloc(sizeof(real_t) * ipm->nb_of_ueq_constr*ipm->optvar_seqlen);
    if (NULL == ipm->P_of_z->P_hat) {return HHMPC_DYNMEM_FAIL;}
    ipm->P_of_z->P_hat_T =
            (real_t *)malloc(sizeof(real_t) * ipm->optvar_seqlen*ipm->nb_of_ueq_constr);
    if (NULL == ipm->P_of_z->P_hat_T) {return HHMPC_DYNMEM_FAIL;}
    ipm->P_of_z->P2_hat =
            (real_t *)malloc(sizeof(real_t) * ipm->nb_of_ueq_constr*ipm->optvar_seqlen);
    if (NULL == ipm->P_of_z->P2_hat) {return HHMPC_DYNMEM_FAIL;}
    ipm->P_of_z->P2_hat_T =
            (real_t *)malloc(sizeof(real_t) * ipm->optvar_seqlen*ipm->nb_of_ueq_constr);
    if (NULL == ipm->P_of_z->P2_hat_T) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->h = ipm->P_of_z->h_hat;
    ipm->P = ipm->P_of_z->P_hat;
    ipm->P_T = ipm->P_of_z->P_hat_T;
    ipm->P2 = ipm->P_of_z->P2_hat;
    ipm->P2_T = ipm->P_of_z->P2_hat_T;
    
    ipm->j_in = &(ipm->conf->in_iter);
    ipm->reg = &(ipm->conf->reg);
    
    ipm->d = (real_t *)malloc(sizeof(real_t) * ipm->nb_of_ueq_constr);
    if (NULL == ipm->d) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->diag_d_sq =
            (real_t *)malloc(sizeof(real_t) * ipm->nb_of_ueq_constr*ipm->nb_of_ueq_constr);
    if (NULL == ipm->diag_d_sq) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->Phi = (real_t *)malloc(sizeof(real_t) * ipm->optvar_seqlen*ipm->optvar_seqlen);
    if (NULL == ipm->Phi) {return HHMPC_DYNMEM_FAIL;}
/*    
    ipm->P_T = (real_t *)malloc(sizeof(real_t) * prb->P->rows*prb->P->cols);
    if (NULL == ipm->P_T) {return HHMPC_DYNMEM_FAIL;}
    mpcinc_mtx_transpose(ipm->P_T, ipm->P, prb->P->rows, prb->P->cols);
*/    
    ipm->A_T = (real_t *)malloc(sizeof(real_t) * prb->A->rows*prb->A->cols);
    if (NULL == ipm->A_T) {return HHMPC_DYNMEM_FAIL;}
    mpcinc_mtx_transpose(ipm->A_T, ipm->A, prb->A->rows, prb->A->cols);
    
    ipm->B_T = (real_t *)malloc(sizeof(real_t) * prb->B->rows*prb->B->cols);
    if (NULL == ipm->B_T) {return HHMPC_DYNMEM_FAIL;}
    mpcinc_mtx_transpose(ipm->B_T, ipm->B, prb->B->rows, prb->B->cols);
    
    ipm->A_B_T =
            (real_t *)malloc(sizeof(real_t) * (prb->A->cols+prb->B->cols)*prb->A->rows);
    if (NULL == ipm->A_B_T) {return HHMPC_DYNMEM_FAIL;}
    
    for (i = 0; i < ipm->state_veclen*ipm->state_veclen; i++)
        ipm->A_B_T[i] = ipm->A_T[i];
    for (i = 0; i < ipm->state_veclen*ipm->control_veclen; i++)
        ipm->A_B_T[ipm->state_veclen*ipm->state_veclen+i] = ipm->B_T[i];
    ipm->A_B =
            (real_t *)malloc(sizeof(real_t) * prb->A->rows*(prb->A->cols+prb->B->cols));
    if (NULL == ipm->A_B) {return HHMPC_DYNMEM_FAIL;}
    mpcinc_mtx_transpose(ipm->A_B, ipm->A_B_T, ipm->optvar_veclen, ipm->state_veclen);
    
    ipm->C_T = (real_t *)malloc(sizeof(real_t) * prb->C->rows*prb->C->cols);
    if (NULL == ipm->C_T) {return HHMPC_DYNMEM_FAIL;}
    mpcinc_mtx_transpose(ipm->C_T, ipm->C, prb->C->rows, prb->C->cols);
    
    ipm->r_d = (real_t *)malloc(ipm->sizeof_optvar_seqlen);
    if (NULL == ipm->r_d) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->r_p = (real_t *)malloc(ipm->sizeof_dual_seqlen);
    if (NULL == ipm->r_p) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->st_size = (real_t *)malloc(sizeof(real_t));
    if (NULL == ipm->st_size) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp1_optvar_seqlen = (real_t *)malloc(ipm->sizeof_optvar_seqlen);
    if (NULL == ipm->tmp1_optvar_seqlen) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp2_optvar_seqlen = (real_t *)malloc(ipm->sizeof_optvar_seqlen);
    if (NULL == ipm->tmp2_optvar_seqlen) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp2_dual_seqlen = (real_t *)malloc(ipm->sizeof_dual_seqlen);
    if (NULL == ipm->tmp2_dual_seqlen) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp3_mtx_optvar_nb_of_ueq = (real_t *)malloc(sizeof(real_t) * ipm->optvar_seqlen*ipm->nb_of_ueq_constr);
    if (NULL == ipm->tmp3_mtx_optvar_nb_of_ueq) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp3_state_veclen = (real_t *)malloc(sizeof(real_t) * ipm->state_veclen);
    if (NULL == ipm->tmp3_state_veclen) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp4_mtx_optvar_optvar = (real_t *)malloc(sizeof(real_t) * ipm->optvar_seqlen*ipm->optvar_seqlen);
    if (NULL == ipm->tmp4_mtx_optvar_optvar) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp4_nb_of_constr = (real_t *)malloc(sizeof(real_t) * ipm->nb_of_ueq_constr);
    if (NULL == ipm->tmp4_nb_of_constr) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp5_nb_of_constr = (real_t *)malloc(sizeof(real_t) * ipm->nb_of_ueq_constr);
    if (NULL == ipm->tmp5_nb_of_constr) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp6_optvar_seqlen = (real_t *)malloc(ipm->sizeof_optvar_seqlen);
    if (NULL == ipm->tmp6_optvar_seqlen) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp7_dual_seqlen = (real_t *)malloc(ipm->sizeof_dual_seqlen);
    if (NULL == ipm->tmp7_dual_seqlen) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp8_L_Y =
            (real_t *)malloc(sizeof(real_t) * (2*ipm->horizon-1)*ipm->state_veclen*ipm->state_veclen);
    if (NULL == ipm->tmp8_L_Y) {return HHMPC_DYNMEM_FAIL;}  /*T*n*n* big enough*/
    
    ipm->tmp9_L_Y_T =
            (real_t *)malloc(sizeof(real_t) * (2*ipm->horizon-1)*ipm->state_veclen*ipm->state_veclen);
    if (NULL == ipm->tmp9_L_Y_T) {return HHMPC_DYNMEM_FAIL;}  /*T*n*n* big enough*/
        
    ipm->tmp8_L_Phi =
            (real_t *)malloc(sizeof(real_t) * ipm->horizon*ipm->optvar_veclen*ipm->optvar_veclen);
    if (NULL == ipm->tmp8_L_Phi) {return HHMPC_DYNMEM_FAIL;}  /*m*m+(T-1)*(n+m)*(n+m)+n*n big enough*/
    
    ipm->tmp9_L_Phi_T =
            (real_t *)malloc(sizeof(real_t) * ipm->horizon*ipm->optvar_veclen*ipm->optvar_veclen);
    if (NULL == ipm->tmp9_L_Phi_T) {return HHMPC_DYNMEM_FAIL;}  /*m*m+(T-1)*(n+m)*(n+m)+n*n big enough*/
    
    ipm->tmp_phibl1 = (real_t *)malloc(sizeof(real_t) * ipm->optvar_veclen*ipm->optvar_veclen);
    if (NULL == ipm->tmp_phibl1) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp_phibl2 = (real_t *)malloc(sizeof(real_t) * ipm->optvar_veclen*ipm->optvar_veclen);
    if (NULL == ipm->tmp_phibl2) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp_phibl3 = (real_t *)malloc(sizeof(real_t) * ipm->optvar_veclen*ipm->optvar_veclen);
    if (NULL == ipm->tmp_phibl3) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp10 = (real_t *)malloc(sizeof(real_t) * ipm->optvar_veclen*ipm->optvar_veclen);
    if (NULL == ipm->tmp10) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmpYbl = (real_t *)malloc(sizeof(real_t) * ipm->state_veclen*ipm->state_veclen);
    if (NULL == ipm->tmpYbl) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmpQbl = (real_t *)malloc(sizeof(real_t) * ipm->state_veclen*ipm->state_veclen);
    if (NULL == ipm->tmpQbl) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp1_res_os = (real_t *)malloc(ipm->sizeof_optvar_seqlen);
    if (NULL == ipm->tmp1_res_os) {return HHMPC_DYNMEM_FAIL;}
    ipm->tmp2_res_os = (real_t *)malloc(ipm->sizeof_optvar_seqlen);
    if (NULL == ipm->tmp2_res_os) {return HHMPC_DYNMEM_FAIL;}
    ipm->tmp3_res_ds = (real_t *)malloc(ipm->sizeof_dual_seqlen);
    if (NULL == ipm->tmp3_res_ds) {return HHMPC_DYNMEM_FAIL;}
        
    ipm->eye_optvar_veclen = (real_t *)malloc(sizeof(real_t) * ipm->optvar_veclen*ipm->optvar_veclen);
    if (NULL == ipm->eye_optvar_veclen) {return HHMPC_DYNMEM_FAIL;}
    eye(ipm->eye_optvar_veclen, ipm->optvar_veclen);
    
    ipm->eye_state_veclen = (real_t *)malloc(sizeof(real_t) * ipm->state_veclen*ipm->state_veclen);
    if (NULL == ipm->eye_state_veclen) {return HHMPC_DYNMEM_FAIL;}
    eye(ipm->eye_state_veclen, ipm->state_veclen);
    
    return HHMPC_DYNMEM_OK;
}

/* Static function definitions */

hhmpc_dynmem_error_t hhmpc_ipm_parse_elements(struct hhmpc_ipm *ipm, cJSON *data)
{
    cJSON *kappa, *optvar, *veclen, *horizon, *state_veclen;
    
    kappa = cJSON_GetObjectItem(data, "kappa");
    if (NULL == kappa) {
        printf("ERROR: could not parse item %s \n", "kappa");
        return HHMPC_DYNMEM_FAIL;
    }
    *(ipm->kappa) = (real_t)kappa->valuedouble;
    
    optvar = cJSON_GetObjectItem(data, "optvar");
    if (NULL == optvar) {
        printf("ERROR: could not parse item %s \n", "optvar");
        return HHMPC_DYNMEM_FAIL;
    }
    veclen = cJSON_GetObjectItem(optvar, "veclen");
    if (NULL == veclen) {
        printf("ERROR: could not parse item %s \n", "veclen");
        return HHMPC_DYNMEM_FAIL;
    }
    ipm->optvar_veclen = (uint32_t)veclen->valueint;
    
    horizon = cJSON_GetObjectItem(optvar, "horizon");
    if (NULL == horizon) {
        printf("ERROR: could not parse item %s \n", "horizon");
        return HHMPC_DYNMEM_FAIL;
    }
    ipm->horizon = (uint32_t)horizon->valueint;
    
    state_veclen = cJSON_GetObjectItem(optvar, "state_veclen");
    if (NULL == state_veclen) {
        printf("ERROR: could not parse item %s \n", "state_veclen");
        return HHMPC_DYNMEM_FAIL;
    }
    ipm->state_veclen = (uint32_t)state_veclen->valueint;
    
    return HHMPC_DYNMEM_OK;
}
