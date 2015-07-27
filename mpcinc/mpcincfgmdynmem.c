#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */

#include "include/cjson.h"
#include "include/mpcincdynmem.h"
#include "include/mpcinccvp.h"
#include "include/mpcincfgmdynmem.h"

/* Static functions declarations */
mpcinc_dynmem_error_t mpcinc_fgm_parse_elements(struct mpcinc_fgm *fgm, cJSON *data);

/* Extern function definitions */

struct mpcinc_fgm *mpcinc_fgm_allocate_solver(void)
{
    struct mpcinc_fgm *fgm = (struct mpcinc_fgm*)malloc(sizeof(struct mpcinc_fgm));
    if (NULL == fgm) {return NULL;}
    struct mpcinc_fgm_conf *conf = (struct mpcinc_fgm_conf*)malloc(sizeof(struct mpcinc_fgm_conf));
    if (NULL == conf) {return NULL;}
    fgm->conf = conf;

    fgm->nu = (real_t *)malloc(sizeof(real_t));
    if (NULL == fgm->nu) {return NULL;}

    fgm->u_ini = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->u_ini) {return NULL;}
    fgm->u_opt = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->u_opt) {return NULL;}

    fgm->tmp1_optvar_seqlen = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->tmp1_optvar_seqlen) {return NULL;}
    fgm->tmp2_optvar_seqlen = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->tmp2_optvar_seqlen) {return NULL;}
    fgm->tmp3_optvar_seqlen = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->tmp3_optvar_seqlen) {return NULL;}
    fgm->tmp4_optvar_seqlen = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->tmp4_optvar_seqlen) {return NULL;}
    fgm->tmp5_optvar_seqlen = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->tmp5_optvar_seqlen) {return NULL;}
    fgm->tmp6_optvar_seqlen = (real_t *)malloc(fgm->sizeof_optvar_seqlen);
    if (NULL == fgm->tmp6_optvar_seqlen) {return NULL;}

    return fgm;
}

mpcinc_dynmem_error_t mpcinc_fgm_setup_solver(struct mpcinc_fgm *fgm, struct mpcinc_cvp_prb *prb, char *fname)
{
    mpcinc_dynmem_error_t ret;
    cJSON *data;

    data = mpcinc_dynmem_get_data(fname);
    if (NULL == data) {return MPCINC_DYNMEM_FAIL;}
    ret = mpcinc_fgm_parse_elements(fgm, data);
    if (MPCINC_DYNMEM_OK != ret) {return ret;}

    fgm->goL = prb->g->data;
    fgm->HoL = prb->H->data;
    fgm->u_lb = prb->u_lb->data;
    fgm->u_ub = prb->u_ub->data;
    fgm->optvar_seqlen = prb->u_lb->rows;
    /* fgm->optvar_veclen = ; FIXME this is needed for warmstart */
    fgm->sizeof_optvar_seqlen = sizeof(real_t) * prb->u_lb->rows;

    fgm->j_in = &(fgm->conf->in_iter);

    return MPCINC_DYNMEM_OK;
}

/* Static function definitions */

mpcinc_dynmem_error_t mpcinc_fgm_parse_elements(struct mpcinc_fgm *fgm, cJSON *data)
{
    cJSON *nu, *optvar, *veclen;

    nu = cJSON_GetObjectItem(data, "nu");
    if (NULL == nu) {
        printf("ERROR: could not parse item %s \n", "nu");
        return MPCINC_DYNMEM_FAIL;
    }
    *(fgm->nu) = (real_t)nu->valuedouble;

    optvar = cJSON_GetObjectItem(data, "optvar");
    if (NULL == optvar) {
        printf("ERROR: could not parse item %s \n", "optvar");
        return MPCINC_DYNMEM_FAIL;
    }
    veclen = cJSON_GetObjectItem(optvar, "veclen");
    if (NULL == veclen) {
        printf("ERROR: could not parse item %s \n", "veclen");
        return MPCINC_DYNMEM_FAIL;
    }
    fgm->optvar_veclen = (uint32_t)veclen->valueint;

    return MPCINC_DYNMEM_OK;
}

