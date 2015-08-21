#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */
#include "include/cjson.h"
#include "include/mpcincdynmem.h"
#include "include/mpcinccvpdynmem.h"

/* Static functions declarations */
static mpcinc_dynmem_error_t mpcinc_get_json_term(struct mpcinc_term *term, cJSON *data, char *jname, char *term_name);
static mpcinc_dynmem_error_t mpcinc_get_json_sub_term(struct mpcinc_term *term, cJSON *data, char *jname, char *term_name, char *sub_name);
static mpcinc_dynmem_error_t mpcinc_parse_elements(struct mpcinc_cvp *cvp, cJSON *data);
static mpcinc_dynmem_error_t mpcinc_get_json_fac_term(struct mpcinc_term *term, cJSON *data, char *jname, char *list_name, int fac_pos);
static mpcinc_dynmem_error_t mpcinc_get_json_term_items(struct mpcinc_term *term, cJSON *jobj);
static mpcinc_dynmem_error_t mpcinc_alloc_data(double **data, int elems);

/* Extern function definitions */

struct mpcinc_cvp *mpcinc_cvp_allocate_former(void) {
    int i;
    struct mpcinc_term *t;
    struct mpcinc_pmetric *p;
    struct mpcinc_cvp *cvp = (struct mpcinc_cvp*)malloc(sizeof(struct mpcinc_cvp));
    if (NULL == cvp) {return NULL;}
    /* parameters */
    t = (struct mpcinc_term*)calloc(MPCINC_PAR_NUM, sizeof(struct mpcinc_term));
    if (NULL == t) {return NULL;}
    for (i=0; i<MPCINC_PAR_NUM; i++) {
        cvp->par[i] = &t[i];
    }

    /* constants */
    t = (struct mpcinc_term*)calloc(MPCINC_CONSTANT_NUM, sizeof(struct mpcinc_term));
    if (NULL == t) {return NULL;}
    for (i=0; i<MPCINC_CONSTANT_NUM; i++) {
        cvp->constant[i] = &t[i];
    }

    /* parametric */
    p = (struct mpcinc_pmetric*)calloc(MPCINC_PMETRIC_NUM, sizeof(struct mpcinc_pmetric));
    if (NULL == p) {return NULL;}
    for (i=0; i<MPCINC_PMETRIC_NUM; i++) {
        cvp->pmetric[i] = &p[i];
        cvp->pmetric[i]->fac_num = (uint32_t*)malloc(sizeof(uint32_t*));
        if (NULL == cvp->pmetric[i]->fac_num) {return NULL;}
        cvp->pmetric[i]->val = (struct mpcinc_term*)malloc(sizeof(struct mpcinc_term));
        if (NULL == cvp->pmetric[i]->val) {return NULL;}
        cvp->pmetric[i]->aux = (struct mpcinc_term*)malloc(sizeof(struct mpcinc_term));
        if (NULL == cvp->pmetric[i]->aux) {return NULL;}
        cvp->pmetric[i]->fac0 = (struct mpcinc_term*)malloc(sizeof(struct mpcinc_term));
        if (NULL == cvp->pmetric[i]->fac0) {return NULL;}
    }

        cvp->pmetric[MPCINC_G]->fac_num[0] = 3;
        cvp->pmetric[MPCINC_G]->fac = (struct mpcinc_term**)calloc(3, sizeof(struct mpcinc_term*));
if (NULL == cvp->pmetric[MPCINC_G]->fac) {return NULL;}
        cvp->pmetric[MPCINC_G]->par = (struct mpcinc_term**)calloc(3, sizeof(struct mpcinc_term*));
if (NULL == cvp->pmetric[MPCINC_G]->par) {return NULL;}
        t = (struct mpcinc_term*)calloc(3, sizeof(struct mpcinc_term));
if (NULL == t) {return NULL;}
        cvp->pmetric[MPCINC_G]->fac[0] = &t[0];
        cvp->pmetric[MPCINC_G]->par[0] = cvp->par[MPCINC_XR];
        cvp->pmetric[MPCINC_G]->fac[1] = &t[1];
        cvp->pmetric[MPCINC_G]->par[1] = cvp->par[MPCINC_X_K];
        cvp->pmetric[MPCINC_G]->fac[2] = &t[2];
        cvp->pmetric[MPCINC_G]->par[2] = cvp->par[MPCINC_UR];


    /* the evaluated problem itself */
    cvp->prb = (struct mpcinc_cvp_prb*)malloc(sizeof(struct mpcinc_cvp_prb));
    if (NULL == cvp->prb) {return NULL;}
        cvp->prb->g = cvp->pmetric[MPCINC_G]->val;
        cvp->prb->H = cvp->constant[MPCINC_H];
        cvp->prb->u_lb = cvp->constant[MPCINC_U_LB];
        cvp->prb->u_ub = cvp->constant[MPCINC_U_UB];


    return cvp;
}

mpcinc_dynmem_error_t mpcinc_cvp_setup_former(struct mpcinc_cvp *cvp, char *fname) {
    mpcinc_dynmem_error_t ret;
    cJSON *data;
    data = mpcinc_dynmem_get_data(fname);
    if (NULL == data) {return MPCINC_DYNMEM_FAIL;}
    ret = mpcinc_parse_elements(cvp, data);
    if (MPCINC_DYNMEM_OK != ret) {return ret;}
    return MPCINC_DYNMEM_OK;
}

/* Static function definitions */
mpcinc_dynmem_error_t mpcinc_parse_elements(struct mpcinc_cvp *cvp, cJSON *data)
{
        mpcinc_get_json_term(cvp->par[MPCINC_XR], data, "par", "xr");
        mpcinc_get_json_term(cvp->par[MPCINC_X_K], data, "par", "x_k");
        mpcinc_get_json_term(cvp->par[MPCINC_UR], data, "par", "ur");
        mpcinc_get_json_sub_term(cvp->pmetric[MPCINC_G]->val, data, "pmetric", "g", "val");
        mpcinc_get_json_sub_term(cvp->pmetric[MPCINC_G]->fac0, data, "pmetric", "g", "fac0");
        mpcinc_get_json_sub_term(cvp->pmetric[MPCINC_G]->aux, data, "pmetric", "g", "aux");
        mpcinc_get_json_fac_term(cvp->pmetric[MPCINC_G]->fac[0], data, "pmetric", "g", 0);
        mpcinc_get_json_fac_term(cvp->pmetric[MPCINC_G]->fac[1], data, "pmetric", "g", 1);
        mpcinc_get_json_fac_term(cvp->pmetric[MPCINC_G]->fac[2], data, "pmetric", "g", 2);
        mpcinc_get_json_term(cvp->constant[MPCINC_H], data, "constant", "H");
        mpcinc_get_json_term(cvp->constant[MPCINC_U_LB], data, "constant", "u_lb");
        mpcinc_get_json_term(cvp->constant[MPCINC_U_UB], data, "constant", "u_ub");

    return MPCINC_DYNMEM_OK;
}

mpcinc_dynmem_error_t mpcinc_get_json_term(struct mpcinc_term *term, cJSON *data, char *jname, char *term_name)
{
    cJSON *jobj = cJSON_GetObjectItem(data, jname);
    if (NULL == jobj) {
        printf("ERROR: could not get item %s \n", jname);
        return MPCINC_DYNMEM_FAIL;
    }
    cJSON *jterm = cJSON_GetObjectItem(jobj, term_name);
    if (NULL == jterm) {
        printf("ERROR: could not get item %s \n", term_name);
        return MPCINC_DYNMEM_FAIL;
    }
    mpcinc_get_json_term_items(term, jterm);
    if (NULL == term) {return MPCINC_DYNMEM_FAIL;}

    return MPCINC_DYNMEM_OK;
 }

mpcinc_dynmem_error_t mpcinc_get_json_sub_term(struct mpcinc_term *term, cJSON *data, char *jname, char *term_name, char *sub_name)
{
    mpcinc_dynmem_error_t ret;

    cJSON *jobj = cJSON_GetObjectItem(data, jname);
    if (NULL == jobj) {
        printf("ERROR: could not get item %s \n", jname);
        return MPCINC_DYNMEM_FAIL;
    }
    cJSON *jpmetric = cJSON_GetObjectItem(jobj, term_name);
    if (NULL == jpmetric) {
        printf("ERROR: could not get item %s \n", term_name);
        return MPCINC_DYNMEM_FAIL;
    }
    cJSON *jterm = cJSON_GetObjectItem(jpmetric, sub_name);
    if (NULL == jterm) {
        printf("ERROR: could not get item %s \n", sub_name);
        return MPCINC_DYNMEM_FAIL;
    }
    ret = mpcinc_get_json_term_items(term, jterm);
    if (MPCINC_DYNMEM_OK != ret) {return ret;}

    return MPCINC_DYNMEM_OK;
 }

mpcinc_dynmem_error_t mpcinc_get_json_fac_term(struct mpcinc_term *term, cJSON *data, char *jname, char *list_name, int fac_pos)
{
    mpcinc_dynmem_error_t ret;

    cJSON *jobj = cJSON_GetObjectItem(data, jname);
    if (NULL == jobj) {
        printf("ERROR: could not get item %s \n", jname);
        return MPCINC_DYNMEM_FAIL;
    }
    cJSON *jpmetric = cJSON_GetObjectItem(jobj, list_name);
    if (NULL == jpmetric) {
        printf("ERROR: could not get item %s \n", list_name);
        return MPCINC_DYNMEM_FAIL;
    }
    cJSON *jfac_list = cJSON_GetObjectItem(jpmetric, "fac");
    if (NULL == jfac_list) {
        printf("ERROR: could not get item %s \n", "fac");
        return MPCINC_DYNMEM_FAIL;
    }
    cJSON *jfac = cJSON_GetArrayItem(jfac_list, fac_pos);
    if (NULL == jfac) {
        printf("ERROR: could not get array item in position %d \n", fac_pos);
        return MPCINC_DYNMEM_FAIL;
    }

    ret = mpcinc_get_json_term_items(term, jfac);
    if (MPCINC_DYNMEM_OK != ret) {return ret;}

    return MPCINC_DYNMEM_OK;
 }

mpcinc_dynmem_error_t mpcinc_get_json_term_items(struct mpcinc_term *term, cJSON *jobj)
{
    cJSON *c;
    mpcinc_dynmem_error_t ret;
    c = cJSON_GetObjectItem(jobj, "cols");
    if (NULL == c) {
        printf("ERROR: could not get item %s \n", "cols");
        return MPCINC_DYNMEM_FAIL;
    }
    term->cols = (uint32_t)c->valueint;
    c = cJSON_GetObjectItem(jobj, "rows");
    if (NULL == c) {
        printf("ERROR: could not get item %s \n", "rows");
        return MPCINC_DYNMEM_FAIL;
    }
    term->rows = (uint32_t)c->valueint;
    cJSON *jdata = cJSON_GetObjectItem(jobj, "data");
    if (NULL == jdata) {
        printf("ERROR: could not get item %s \n", "data");
        return MPCINC_DYNMEM_FAIL;
    }
    int elems = (int)(term->rows * term->cols);
    if (elems != cJSON_GetArraySize(jdata)) {
        printf("Size error, expected: %d*%d (rows*cols); actual: %d \n",
            term->rows, term->cols, cJSON_GetArraySize(jdata));
        {return MPCINC_DYNMEM_FAIL;}
    } else {
        ret = mpcinc_alloc_data(&(term->data), elems);
        if (MPCINC_DYNMEM_OK != ret) {return ret;}
        int i;
        for (i=0;i<elems;i++) {
            term->data[i] = cJSON_GetArrayItem(jdata, i)->valuedouble;
        }
    }

    return MPCINC_DYNMEM_OK;
 }

mpcinc_dynmem_error_t mpcinc_alloc_data(double **data, int elems)
{
    data[0] = (double*) malloc(elems*sizeof(double));
    if (NULL == data[0]) {return MPCINC_DYNMEM_FAIL;}

    return MPCINC_DYNMEM_OK;
}
