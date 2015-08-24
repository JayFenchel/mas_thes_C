#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */

#include "mc04types.h"
#include "include/cjson.h"
#include "include/hhmpcipmdynmem.h"

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
    
    data = hhmpc_dynmem_get_data(fname);
    if (NULL == data) {return HHMPC_DYNMEM_FAIL;}
    ret = hhmpc_ipm_parse_elements(ipm, data);
    if (HHMPC_DYNMEM_OK != ret) {return ret;}
    
    ipm->b = prb->b->data;
    ipm->h = prb->h->data;
    ipm->q = prb->q->data;
    ipm->r = prb->r->data;
    ipm->C = prb->C->data;
    
    ipm->optvar_seqlen = ipm->optvar_veclen * ipm->horizon;
    ipm->optvar_dual = ipm->dim_state * ipm->horizon;
    ipm->sizeof_optvar_seqlen = sizeof(real_t) * ipm->optvar_seqlen;
    ipm->sizeof_optvar_dual = sizeof(real_t) * ipm->optvar_dual;
    
    ipm->j_in = &(ipm->conf->in_iter);
    
    ipm->z_ini = (real_t *)malloc(ipm->sizeof_optvar_seqlen);
    if (NULL == ipm->z_ini) {return HHMPC_DYNMEM_FAIL;}
    ipm->z_opt = (real_t *)malloc(ipm->sizeof_optvar_seqlen);
    if (NULL == ipm->z_opt) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->r_p = (real_t *)malloc(ipm->sizeof_optvar_dual);
    if (NULL == ipm->r_p) {return HHMPC_DYNMEM_FAIL;}
    
    ipm->tmp1_optvar_seqlen = (real_t *)malloc(ipm->sizeof_optvar_seqlen);
    if (NULL == ipm->tmp1_optvar_seqlen) {return HHMPC_DYNMEM_FAIL;}

    return HHMPC_DYNMEM_OK;
}

/* Static function definitions */

hhmpc_dynmem_error_t hhmpc_ipm_parse_elements(struct hhmpc_ipm *ipm, cJSON *data)
{
    cJSON *kappa, *optvar, *veclen, *horizon, *dim_state;
    
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
    
    dim_state = cJSON_GetObjectItem(optvar, "dim_state");
    if (NULL == dim_state) {
        printf("ERROR: could not parse item %s \n", "dim_state");
        return HHMPC_DYNMEM_FAIL;
    }
    ipm->dim_state = (uint32_t)dim_state->valueint;
    
    return HHMPC_DYNMEM_OK;
}
