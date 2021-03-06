#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */
#include "include/cjson.h"

#include "include/hhmpcsocpdynmem.h"
#include "include/mpcincmtxops.h"

/* Static functions declarations */
static hhmpc_dynmem_error_t hhmpc_get_json_term(struct hhmpc_term *term,
                                                cJSON *data,
                                                char *jname,
                                                char *term_name);
static hhmpc_dynmem_error_t hhmpc_get_json_sub_term(struct hhmpc_term *term,
                                                    cJSON *data,
                                                    char *jname,
                                                    char *term_name,
                                                    char *sub_name);
static hhmpc_dynmem_error_t hhmpc_get_json_constr_term(struct hhmpc_term *term,
                                                       cJSON *data,
                                                       char *jname,
                                                       char *list_name,
                                                       int fac_pos);
static hhmpc_dynmem_error_t hhmpc_get_json_fac_term(struct hhmpc_term *term,
                                                    cJSON *data,
                                                    char *jname,
                                                    char *list_name,
                                                    int fac_pos);
static hhmpc_dynmem_error_t hhmpc_parse_elements(struct hhmpc_socp *socp,
                                                 cJSON *data);
static hhmpc_dynmem_error_t hhmpc_get_json_term_items(struct hhmpc_term *term,
                                                      cJSON *jobj);
static hhmpc_dynmem_error_t hhmpc_alloc_data(double **data, int elems);

/* Extern function definitions */

struct hhmpc_socp *hhmpc_socp_allocate_former(void)
{
    int i;
    struct hhmpc_term *t;
    struct hhmpc_pmetric *p;
    struct hhmpc_socp *socp = (struct hhmpc_socp*)malloc(sizeof(struct hhmpc_socp));
    if (NULL == socp) {return NULL;}
    
    /* parameters */
    t = (struct hhmpc_term*)calloc(HHMPC_PAR_NUM, sizeof(struct hhmpc_term));
    if (NULL == t) {return NULL;}
    for (i=0; i<HHMPC_PAR_NUM; i++) {
        socp->par[i] = &t[i];
    }
    
    /* constants */
    t = (struct hhmpc_term*)calloc(HHMPC_CONST_NUM, sizeof(struct hhmpc_term));
    if (NULL == t) {return NULL;}
    for (i=0; i<HHMPC_CONST_NUM; i++) {
        socp->constant[i] = &t[i];
    }
    
    /* parametric */
    p = (struct hhmpc_pmetric*)calloc(HHMPC_PMETRIC_NUM, sizeof(struct hhmpc_pmetric));
    if (NULL == p) {return NULL;}
    for (i=0; i<HHMPC_PMETRIC_NUM; i++) {
        socp->pmetric[i] = &p[i];
        socp->pmetric[i]->fac_num = (uint32_t*)malloc(sizeof(uint32_t)); /* statt sizeof(uint32_t*) */
        if (NULL == socp->pmetric[i]->fac_num) {return NULL;}
        socp->pmetric[i]->val = (struct hhmpc_term*)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->pmetric[i]->val) {return NULL;}
        socp->pmetric[i]->aux = (struct hhmpc_term*)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->pmetric[i]->aux) {return NULL;}
        socp->pmetric[i]->fac0 = (struct hhmpc_term*)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->pmetric[i]->fac0) {return NULL;}
    }
    
    /* pmetric ZR allocieren und mit richtigem par verknüpfen*/
    socp->pmetric[HHMPC_ZR]->fac_num[0] = 2;
    socp->pmetric[HHMPC_ZR]->fac =
            (struct hhmpc_term**)calloc(2, sizeof(struct hhmpc_term*));
    if (NULL == socp->pmetric[HHMPC_ZR]->fac) {return NULL;}
    socp->pmetric[HHMPC_ZR]->par =
            (struct hhmpc_term**)calloc(2, sizeof(struct hhmpc_term*));
    if (NULL == socp->pmetric[HHMPC_ZR]->par) {return NULL;}
    t = (struct hhmpc_term*)calloc(2, sizeof(struct hhmpc_term));
    if (NULL == t) {return NULL;}
    
    socp->pmetric[HHMPC_ZR]->fac[0] = &t[0];
    socp->pmetric[HHMPC_ZR]->par[0] = socp->par[HHMPC_UR];
    socp->pmetric[HHMPC_ZR]->fac[1] = &t[1];
    socp->pmetric[HHMPC_ZR]->par[1] = socp->par[HHMPC_XR];
    
    /* pmetric B_KL allocieren und mit richtigem par verknüpfen*/
    socp->pmetric[HHMPC_B_KL]->fac_num[0] = 1;
    socp->pmetric[HHMPC_B_KL]->fac =
            (struct hhmpc_term**)calloc(1, sizeof(struct hhmpc_term*));
    if (NULL == socp->pmetric[HHMPC_B_KL]->fac) {return NULL;}
    socp->pmetric[HHMPC_B_KL]->par =
            (struct hhmpc_term**)calloc(1, sizeof(struct hhmpc_term*));
    if (NULL == socp->pmetric[HHMPC_B_KL]->par) {return NULL;}
    t = (struct hhmpc_term*)calloc(1, sizeof(struct hhmpc_term));
    if (NULL == t) {return NULL;}
    
    socp->pmetric[HHMPC_B_KL]->fac[0] = &t[0];
    socp->pmetric[HHMPC_B_KL]->par[0] = socp->par[HHMPC_XK];
    
    /* pmetric H_KL allocieren und mit richtigem par verknüpfen*/
    socp->pmetric[HHMPC_H_KL]->fac_num[0] = 1;
    socp->pmetric[HHMPC_H_KL]->fac =
            (struct hhmpc_term**)calloc(1, sizeof(struct hhmpc_term*));
    if (NULL == socp->pmetric[HHMPC_H_KL]->fac) {return NULL;}
    socp->pmetric[HHMPC_H_KL]->par =
            (struct hhmpc_term**)calloc(1, sizeof(struct hhmpc_term*));
    if (NULL == socp->pmetric[HHMPC_H_KL]->par) {return NULL;}
    t = (struct hhmpc_term*)calloc(1, sizeof(struct hhmpc_term));
    if (NULL == t) {return NULL;}
    
    socp->pmetric[HHMPC_H_KL]->fac[0] = &t[0];
    socp->pmetric[HHMPC_H_KL]->par[0] = socp->par[HHMPC_XK];
    
    /* pmetric HSOFT_KL allocieren und mit richtigem par verknüpfen*/
    socp->pmetric[HHMPC_HSOFT_KL]->fac_num[0] = 1;
    socp->pmetric[HHMPC_HSOFT_KL]->fac =
            (struct hhmpc_term**)calloc(1, sizeof(struct hhmpc_term*));
    if (NULL == socp->pmetric[HHMPC_HSOFT_KL]->fac) {return NULL;}
    socp->pmetric[HHMPC_HSOFT_KL]->par =
            (struct hhmpc_term**)calloc(1, sizeof(struct hhmpc_term*));
    if (NULL == socp->pmetric[HHMPC_HSOFT_KL]->par) {return NULL;}
    t = (struct hhmpc_term*)calloc(1, sizeof(struct hhmpc_term));
    if (NULL == t) {return NULL;}
    
    socp->pmetric[HHMPC_HSOFT_KL]->fac[0] = &t[0];
    socp->pmetric[HHMPC_HSOFT_KL]->par[0] = socp->par[HHMPC_XK];
    
    /* pmetric G_KL allocieren und mit richtigem par verknüpfen*/
    socp->pmetric[HHMPC_G_KL]->fac_num[0] = 1;
    socp->pmetric[HHMPC_G_KL]->fac =
            (struct hhmpc_term**)calloc(1, sizeof(struct hhmpc_term*));
    if (NULL == socp->pmetric[HHMPC_G_KL]->fac) {return NULL;}
    socp->pmetric[HHMPC_G_KL]->par =
            (struct hhmpc_term**)calloc(1, sizeof(struct hhmpc_term*));
    if (NULL == socp->pmetric[HHMPC_G_KL]->par) {return NULL;}
    t = (struct hhmpc_term*)calloc(1, sizeof(struct hhmpc_term));
    if (NULL == t) {return NULL;}
    
    socp->pmetric[HHMPC_G_KL]->fac[0] = &t[0];
    socp->pmetric[HHMPC_G_KL]->par[0] = socp->par[HHMPC_XK];
    
    
    /* the evaluated problem itself */
    socp->prb = (struct hhmpc_socp_prb*)malloc(sizeof(struct hhmpc_socp_prb));
    if (NULL == socp->prb) {return NULL;}
        socp->prb->z_ini = socp->par[HHMPC_ZINI];
        socp->prb->v_ini = socp->par[HHMPC_VINI];
        socp->prb->x_k = socp->par[HHMPC_XK];
        socp->prb->tmp_state_veclen = 
                (struct hhmpc_term *)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->prb->tmp_state_veclen) {return NULL;}
        socp->prb->u_k = 
                (struct hhmpc_term *)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->prb->u_k) {return NULL;}
        socp->prb->Psoft =
                (struct hhmpc_term *)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->prb->Psoft) {return NULL;}
    
        socp->prb->zref = socp->pmetric[HHMPC_ZR]->val;
        socp->prb->b = socp->pmetric[HHMPC_B_KL]->val;
        socp->prb->h = socp->pmetric[HHMPC_H_KL]->val;
        socp->prb->hsoft = socp->pmetric[HHMPC_HSOFT_KL]->val;
        socp->prb->Fusoft = socp->constant[HHMPC_FUSOFT];
        socp->prb->Fxsoft = socp->constant[HHMPC_FXSOFT];
        socp->prb->Ffsoft = socp->constant[HHMPC_FFSOFT];
        socp->prb->Psoft = socp->constant[HHMPC_PSOFT];        
        socp->prb->g = socp->pmetric[HHMPC_G_KL]->val;
        socp->prb->q = socp->constant[HHMPC_Q_KL];
        socp->prb->r = socp->constant[HHMPC_R_KL];
        socp->prb->S = socp->constant[HHMPC_S];
        socp->prb->S_T = socp->constant[HHMPC_S_T];
        socp->prb->A = socp->constant[HHMPC_A];
        socp->prb->B = socp->constant[HHMPC_B];
        socp->prb->C = socp->constant[HHMPC_C];
        socp->prb->P = socp->constant[HHMPC_P];
        socp->prb->H = socp->constant[HHMPC_H];
    
    return socp;
}

hhmpc_dynmem_error_t hhmpc_socp_setup_former(struct hhmpc_socp *socp,
                                             char *fname)
{
    hhmpc_dynmem_error_t ret;
    cJSON *data;
    data = hhmpc_dynmem_get_data(fname);
    if (NULL == data) {return HHMPC_DYNMEM_FAIL;}
    ret = hhmpc_parse_elements(socp, data);    
    if(HHMPC_DYNMEM_OK != ret) {return ret;}
    
    return HHMPC_DYNMEM_OK;
}

/* Static function definitions */

hhmpc_dynmem_error_t hhmpc_parse_elements(struct hhmpc_socp *socp, cJSON *data)
{
    struct hhmpc_socc *s;
    struct hhmpc_qc *q;
    uint32_t i, j;
    /* get optvar veclen an horizon to calculate length of g and ... */
    cJSON *optvar, *veclen, *horizon, *nb_socc, *nb_qc, *jobj, *jconstr, *c;
    
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
    socp->prb->optvar_veclen = (uint32_t)veclen->valueint;
    
    horizon = cJSON_GetObjectItem(optvar, "horizon");
    if (NULL == horizon) {
        printf("ERROR: could not parse item %s \n", "horizon");
        return HHMPC_DYNMEM_FAIL;
    }
    socp->prb->horizon = (uint32_t)horizon->valueint;
    
    socp->prb->optvar_seqlen = socp->prb->optvar_veclen * socp->prb->horizon;
    socp->prb->sizeof_optvar_seqlen = sizeof(real_t) * socp->prb->optvar_seqlen;
    
    nb_socc = cJSON_GetObjectItem(data, "nb_socc");
    if (NULL == nb_socc) {
        printf("ERROR: could not parse item %s \n", "nb_socc");
        return HHMPC_DYNMEM_FAIL;
    }
    socp->prb->nb_socc = (uint32_t)nb_socc->valueint;
    socp->prb->socc = 
            (struct hhmpc_socc**)calloc(socp->prb->nb_socc, sizeof(struct hhmpc_socc*));
    if (NULL == socp->prb->socc) {return HHMPC_DYNMEM_FAIL;}
    s = (struct hhmpc_socc*)calloc(socp->prb->nb_socc, sizeof(struct hhmpc_socc));
    if (NULL == s) {return HHMPC_DYNMEM_FAIL;}
    for (i = 0; i < socp->prb->nb_socc; i++){
        socp->prb->socc[i] = &s[i];
        socp->prb->socc[i]->A = (struct hhmpc_term*)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->prb->socc[i]->A) {return HHMPC_DYNMEM_FAIL;}
        hhmpc_get_json_constr_term(socp->prb->socc[i]->A, data, "socc", "A", i);
        socp->prb->socc[i]->b = (struct hhmpc_term*)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->prb->socc[i]->b) {return HHMPC_DYNMEM_FAIL;}
        hhmpc_get_json_constr_term(socp->prb->socc[i]->b, data, "socc", "b", i);
        socp->prb->socc[i]->c = (struct hhmpc_term*)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->prb->socc[i]->c) {return HHMPC_DYNMEM_FAIL;}
        hhmpc_get_json_constr_term(socp->prb->socc[i]->c, data, "socc", "c", i);
        socp->prb->socc[i]->d = (struct hhmpc_term*)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->prb->socc[i]->d) {return HHMPC_DYNMEM_FAIL;}
        hhmpc_get_json_constr_term(socp->prb->socc[i]->d, data, "socc", "d", i);
        jobj = cJSON_GetObjectItem(data, "socc");
        if (NULL == jobj) {
            printf("ERROR: could not get item %s \n", "socc");
            return HHMPC_DYNMEM_FAIL;
        }
        jconstr = cJSON_GetArrayItem(jobj, i);
        if (NULL == jconstr) {
            printf("ERROR: could not get array item in position %d \n", i);
            return HHMPC_DYNMEM_FAIL;
        }
        c = cJSON_GetObjectItem(jconstr, "par0");
        if (NULL == c) {
            printf("ERROR: could not get item %s \n", "par0");
            return HHMPC_DYNMEM_FAIL;
        }
        socp->prb->socc[i]->par_0 = (uint32_t)c->valueint;
        c = cJSON_GetObjectItem(jconstr, "parl");
        if (NULL == c) {
            printf("ERROR: could not get item %s \n", "parl");
            return HHMPC_DYNMEM_FAIL;
        }
        socp->prb->socc[i]->par_l = (uint32_t)c->valueint;
    }
    
    nb_qc = cJSON_GetObjectItem(data, "nb_qc");
    if (NULL == nb_qc) {
        printf("ERROR: could not parse item %s \n", "nb_qc");
        return HHMPC_DYNMEM_FAIL;
    }
    socp->prb->nb_qc = (uint32_t)nb_qc->valueint;
    socp->prb->qc = 
            (struct hhmpc_qc**)calloc(socp->prb->nb_qc, sizeof(struct hhmpc_qc*));
    if (NULL == socp->prb->qc) {return HHMPC_DYNMEM_FAIL;}
    q = (struct hhmpc_qc*)calloc(socp->prb->nb_qc, sizeof(struct hhmpc_qc));
    if (NULL == q) {return HHMPC_DYNMEM_FAIL;}
    for (i = 0; i < socp->prb->nb_qc; i++){
        socp->prb->qc[i] = &q[i];
        socp->prb->qc[i]->Gamma = (struct hhmpc_term*)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->prb->qc[i]->Gamma) {return HHMPC_DYNMEM_FAIL;}
        hhmpc_get_json_constr_term(socp->prb->qc[i]->Gamma, data, "qc", "Gamma", i);
        socp->prb->qc[i]->beta = (struct hhmpc_term*)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->prb->qc[i]->beta) {return HHMPC_DYNMEM_FAIL;}
        hhmpc_get_json_constr_term(socp->prb->qc[i]->beta , data, "qc", "beta", i);
        socp->prb->qc[i]->alpha = (struct hhmpc_term*)malloc(sizeof(struct hhmpc_term));
        if (NULL == socp->prb->qc[i]->alpha) {return HHMPC_DYNMEM_FAIL;}
        hhmpc_get_json_constr_term(socp->prb->qc[i]->alpha, data, "qc", "alpha", i);
        jobj = cJSON_GetObjectItem(data, "qc");
        if (NULL == jobj) {
            printf("ERROR: could not get item %s \n", "qc");
            return HHMPC_DYNMEM_FAIL;
        }
        jconstr = cJSON_GetArrayItem(jobj, i);
        if (NULL == jconstr) {
            printf("ERROR: could not get array item in position %d \n", i);
            return HHMPC_DYNMEM_FAIL;
        }
        c = cJSON_GetObjectItem(jconstr, "par0");
        if (NULL == c) {
            printf("ERROR: could not get item %s \n", "par0");
            return HHMPC_DYNMEM_FAIL;
        }
        socp->prb->qc[i]->par_0 = (uint32_t)c->valueint;
        c = cJSON_GetObjectItem(jconstr, "parl");
        if (NULL == c) {
            printf("ERROR: could not get item %s \n", "parl");
            return HHMPC_DYNMEM_FAIL;
        }
        socp->prb->qc[i]->par_l = (uint32_t)c->valueint;
    }
    
    hhmpc_get_json_term(socp->par[HHMPC_ZINI], data, "par", "zini");
    hhmpc_get_json_term(socp->par[HHMPC_VINI], data, "par", "vini");
    hhmpc_get_json_term(socp->par[HHMPC_UR], data, "par", "ur");
    hhmpc_get_json_term(socp->par[HHMPC_XR], data, "par", "xr");
    hhmpc_get_json_term(socp->par[HHMPC_XK], data, "par", "xk");
    socp->prb->tmp_state_veclen->data =
                (real_t *)malloc(sizeof(real_t) * socp->prb->x_k->rows);
    if (NULL == socp->prb->tmp_state_veclen->data) {return HHMPC_DYNMEM_FAIL;}
                
    socp->prb->tmp_state_veclen->rows = socp->prb->x_k->rows;
    socp->prb->tmp_state_veclen->cols = socp->prb->x_k->cols;
    
    hhmpc_get_json_sub_term(socp->pmetric[HHMPC_ZR]->val, data, "pmetric", "zref", "val");
    hhmpc_get_json_sub_term(socp->pmetric[HHMPC_ZR]->fac0, data, "pmetric", "zref", "fac0");
    hhmpc_get_json_sub_term(socp->pmetric[HHMPC_ZR]->aux, data, "pmetric", "zref", "aux");
    hhmpc_get_json_fac_term(socp->pmetric[HHMPC_ZR]->fac[0], data, "pmetric", "zref", 0);
    hhmpc_get_json_fac_term(socp->pmetric[HHMPC_ZR]->fac[1], data, "pmetric", "zref", 1);
    
    hhmpc_get_json_sub_term(socp->pmetric[HHMPC_B_KL]->val, data, "pmetric", "b", "val");
    hhmpc_get_json_sub_term(socp->pmetric[HHMPC_B_KL]->fac0, data, "pmetric", "b", "fac0");
    hhmpc_get_json_sub_term(socp->pmetric[HHMPC_B_KL]->aux, data, "pmetric", "b", "aux");
    hhmpc_get_json_fac_term(socp->pmetric[HHMPC_B_KL]->fac[0], data, "pmetric", "b", 0);
    
    hhmpc_get_json_sub_term(socp->pmetric[HHMPC_H_KL]->val, data, "pmetric", "h", "val");
    hhmpc_get_json_sub_term(socp->pmetric[HHMPC_H_KL]->fac0, data, "pmetric", "h", "fac0");
    hhmpc_get_json_sub_term(socp->pmetric[HHMPC_H_KL]->aux, data, "pmetric", "h", "aux");
    hhmpc_get_json_fac_term(socp->pmetric[HHMPC_H_KL]->fac[0], data, "pmetric", "h", 0);
    
    hhmpc_get_json_sub_term(socp->pmetric[HHMPC_HSOFT_KL]->val, data, "pmetric", "hsoft", "val");
    hhmpc_get_json_sub_term(socp->pmetric[HHMPC_HSOFT_KL]->fac0, data, "pmetric", "hsoft", "fac0");
    hhmpc_get_json_sub_term(socp->pmetric[HHMPC_HSOFT_KL]->aux, data, "pmetric", "hsoft", "aux");
    hhmpc_get_json_fac_term(socp->pmetric[HHMPC_HSOFT_KL]->fac[0], data, "pmetric", "hsoft", 0);
    
    hhmpc_get_json_term(socp->constant[HHMPC_Q_KL], data, "constant", "q");
    hhmpc_get_json_term(socp->constant[HHMPC_R_KL], data, "constant", "r");
    hhmpc_get_json_term(socp->constant[HHMPC_S], data, "constant", "S");
    hhmpc_get_json_term(socp->constant[HHMPC_A], data, "constant", "A");
    hhmpc_get_json_term(socp->constant[HHMPC_B], data, "constant", "B");
    hhmpc_get_json_term(socp->constant[HHMPC_C], data, "constant", "C");
    hhmpc_get_json_term(socp->constant[HHMPC_P], data, "constant", "P");
    hhmpc_get_json_term(socp->constant[HHMPC_H], data, "constant", "H");
    
    /* Form Psoft_T and hsoft*/
    hhmpc_get_json_term(socp->constant[HHMPC_FUSOFT], data, "constant", "Fusoft");
    hhmpc_get_json_term(socp->constant[HHMPC_FXSOFT], data, "constant", "Fxsoft");
    hhmpc_get_json_term(socp->constant[HHMPC_FFSOFT], data, "constant", "Ffsoft");
    hhmpc_get_json_term(socp->constant[HHMPC_PSOFT], data, "constant", "Psoft");
//     struct hhmpc_term *Ps = socp->prb->Psoft;
//     Ps->rows = socp->prb->horizon*socp->constant[HHMPC_FUSOFT]->rows +
//             socp->constant[HHMPC_FFSOFT]->rows;
//     Ps->cols = socp->prb->optvar_seqlen;
    
 /* Values of pmetric HHMPC_G_KL are set here */
    socp->prb->S_T->rows = socp->constant[HHMPC_S]->cols;
    socp->prb->S_T->cols = socp->constant[HHMPC_S]->rows;
    socp->prb->S_T->data =
            (real_t *)malloc(sizeof(real_t) * socp->prb->S_T->rows * socp->prb->S_T->cols);
    if (NULL == socp->prb->S_T->data) {return HHMPC_DYNMEM_FAIL;}
    mpcinc_mtx_transpose(socp->prb->S_T->data, socp->constant[HHMPC_S]->data,
                         socp->constant[HHMPC_S]->rows, socp->constant[HHMPC_S]->cols);
 
    socp->pmetric[HHMPC_G_KL]->val->data =
            (real_t *)malloc(socp->prb->sizeof_optvar_seqlen);
    if (NULL == socp->pmetric[HHMPC_G_KL]->val->data) {return HHMPC_DYNMEM_FAIL;}
    socp->pmetric[HHMPC_G_KL]->val->rows = socp->prb->optvar_seqlen;
    socp->pmetric[HHMPC_G_KL]->val->cols = 1;
    socp->pmetric[HHMPC_G_KL]->aux->data =
            (real_t *)malloc(sizeof(real_t) * socp->constant[HHMPC_R_KL]->rows);
    if (NULL == socp->pmetric[HHMPC_G_KL]->aux->data) {return HHMPC_DYNMEM_FAIL;}
    socp->pmetric[HHMPC_G_KL]->aux->rows = socp->constant[HHMPC_R_KL]->rows;
    socp->pmetric[HHMPC_G_KL]->aux->cols = 1;
    socp->pmetric[HHMPC_G_KL]->fac[0]->data =
            (real_t *)malloc(sizeof(real_t) * socp->prb->S_T->rows * socp->prb->S_T->cols);
    if (NULL == socp->pmetric[HHMPC_G_KL]->fac[0]->data) {return HHMPC_DYNMEM_FAIL;}
    mpcinc_mtx_scale(socp->pmetric[HHMPC_G_KL]->fac[0]->data, socp->prb->S_T->data, 2.,
                     socp->prb->S_T->rows, socp->prb->S_T->cols);
    socp->pmetric[HHMPC_G_KL]->fac[0]->rows = socp->prb->S_T->rows;
    socp->pmetric[HHMPC_G_KL]->fac[0]->cols = socp->prb->S_T->cols;
    
    socp->pmetric[HHMPC_G_KL]->fac0->data =
            (real_t *)malloc(socp->prb->sizeof_optvar_seqlen);
    for (i = 0; i < socp->prb->horizon; i++){
        for (j = 0; j < socp->constant[HHMPC_R_KL]->rows; j++){
            socp->pmetric[HHMPC_G_KL]->fac0->data[i*socp->prb->optvar_veclen + j] =
                    socp->constant[HHMPC_R_KL]->data[j];
        }
        for (j = 0; j < socp->constant[HHMPC_Q_KL]->rows; j++){
            socp->pmetric[HHMPC_G_KL]->fac0->data[i*socp->prb->optvar_veclen + socp->constant[HHMPC_R_KL]->rows + j] =
                    socp->constant[HHMPC_Q_KL]->data[j];
        }
    }
    
    return HHMPC_DYNMEM_OK;
}

hhmpc_dynmem_error_t hhmpc_get_json_term (struct hhmpc_term *term, cJSON *data,
                                          char *jname, char *term_name )
{
    cJSON *jobj = cJSON_GetObjectItem(data, jname);
    if (NULL == jobj) {
        printf("ERROR: could not get item %s \n", jname);
        return HHMPC_DYNMEM_FAIL;
    }
    cJSON *jterm = cJSON_GetObjectItem(jobj, term_name);
    if (NULL == jterm) {
        printf("ERROR: could not get item %s \n", term_name);
        return HHMPC_DYNMEM_FAIL;
    }
    hhmpc_get_json_term_items(term, jterm);
    if (NULL == term) {return HHMPC_DYNMEM_FAIL;}

    return HHMPC_DYNMEM_OK;
}

hhmpc_dynmem_error_t hhmpc_get_json_sub_term(struct hhmpc_term *term,
                                             cJSON *data, char *jname,
                                             char *term_name, char *sub_name)
{
    hhmpc_dynmem_error_t ret;
    cJSON *jobj = cJSON_GetObjectItem(data, jname);
    if (NULL == jobj) {
        printf("ERROR: could not get item %s \n", jname);
        return HHMPC_DYNMEM_FAIL;
    }
    cJSON *jpmetric = cJSON_GetObjectItem(jobj, term_name);
    if (NULL == jpmetric) {
        printf("ERROR: could not get item %s \n", term_name);
        return HHMPC_DYNMEM_FAIL;
    }
    cJSON *jterm = cJSON_GetObjectItem(jpmetric, sub_name);
    if (NULL == jterm) {
        printf("ERROR: could not get item %s \n", sub_name);
        return HHMPC_DYNMEM_FAIL;
    }
    ret = hhmpc_get_json_term_items(term, jterm);
    if (HHMPC_DYNMEM_OK != ret) {return ret;}
    
    return HHMPC_DYNMEM_OK;
}

hhmpc_dynmem_error_t hhmpc_get_json_constr_term(struct hhmpc_term *term,
                                                cJSON *data, char *jname,
                                                char *list_name, int fac_pos)
{
    hhmpc_dynmem_error_t ret;
    cJSON *jobj = cJSON_GetObjectItem(data, jname);
    if (NULL == jobj) {
        printf("ERROR: could not get item %s \n", jname);
        return HHMPC_DYNMEM_FAIL;
    }
    cJSON *jconstr = cJSON_GetArrayItem(jobj, fac_pos);
    if (NULL == jconstr) {
        printf("ERROR: could not get array item in position %d \n", fac_pos);
        return HHMPC_DYNMEM_FAIL;
    }
    cJSON *jc_mtx = cJSON_GetObjectItem(jconstr, list_name);
    if (NULL == jc_mtx) {
        printf("ERROR: could not get item %s \n", list_name);
        return HHMPC_DYNMEM_FAIL;
    }
    
    ret = hhmpc_get_json_term_items(term, jc_mtx);
    if (HHMPC_DYNMEM_OK != ret) {return ret;}
    
    return HHMPC_DYNMEM_OK;
}

hhmpc_dynmem_error_t hhmpc_get_json_fac_term(struct hhmpc_term *term,
                                             cJSON *data, char *jname,
                                             char *list_name, int fac_pos)
{
    hhmpc_dynmem_error_t ret;
    cJSON *jobj = cJSON_GetObjectItem(data, jname);
    if (NULL == jobj) {
        printf("ERROR: could not get item %s \n", jname);
        return HHMPC_DYNMEM_FAIL;
    }
    cJSON *jpmetric = cJSON_GetObjectItem(jobj, list_name);
    if (NULL == jpmetric) {
        printf("ERROR: could not get item %s \n", list_name);
        return HHMPC_DYNMEM_FAIL;
    }
    cJSON *jfac_list = cJSON_GetObjectItem(jpmetric, "fac");
    if (NULL == jfac_list) {
        printf("ERROR: could not get item %s \n", "fac");
        return HHMPC_DYNMEM_FAIL;
    }
    cJSON *jfac = cJSON_GetArrayItem(jfac_list, fac_pos);
    if (NULL == jfac) {
        printf("ERROR: could not get array item in position %d \n", fac_pos);
        return HHMPC_DYNMEM_FAIL;
    }
    
    ret = hhmpc_get_json_term_items(term, jfac);
    if (HHMPC_DYNMEM_OK != ret) {return ret;}

    return HHMPC_DYNMEM_OK;
}


hhmpc_dynmem_error_t hhmpc_get_json_term_items(struct hhmpc_term *term,
                                               cJSON *jobj)
{
    cJSON *c;
    hhmpc_dynmem_error_t ret;
    c = cJSON_GetObjectItem(jobj, "cols");
    if (NULL == c) {
        printf("ERROR: could not get item %s \n", "cols");
        return HHMPC_DYNMEM_FAIL;
    }
    term->cols = (uint32_t)c->valueint;
    c = cJSON_GetObjectItem(jobj, "rows");
    if (NULL == c) {
        printf("ERROR: could not get item %s \n", "rows");
        return HHMPC_DYNMEM_FAIL;
    }
    term->rows = (uint32_t)c->valueint;
    cJSON *jdata = cJSON_GetObjectItem(jobj, "data");
    if (NULL == jdata) {
        printf("ERROR: could not get item %s \n", "data");
        return HHMPC_DYNMEM_FAIL;
    }
    int elems = (int)(term->rows * term->cols);
    if (elems != cJSON_GetArraySize(jdata)) {
        printf("Size error, expected: %d*%d (rows*cols); actual: %d \n",
            term->rows, term->cols, cJSON_GetArraySize(jdata));
        {return HHMPC_DYNMEM_FAIL;}
    } else {
        ret = hhmpc_alloc_data(&(term->data), elems);
        if (HHMPC_DYNMEM_OK != ret) {return ret;}
        int i;
        for (i=0;i<elems;i++) {
            term->data[i] = cJSON_GetArrayItem(jdata, i)->valuedouble;
        }
    }
    
    return HHMPC_DYNMEM_OK;
}

hhmpc_dynmem_error_t hhmpc_alloc_data (double **data, int elems)
{
    data[0] = (double*) malloc(elems*sizeof(double));
    if (NULL == data[0]) {return HHMPC_DYNMEM_FAIL;}
    
    return HHMPC_DYNMEM_OK;
}
