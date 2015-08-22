#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */
#include "include/cjson.h"

#include "include/hhmpcsocpdynmem.h"

/* Static functions declarations */
static hhmpc_dynmem_error_t hhmpc_get_json_term(struct hhmpc_term *term,
                                                cJSON *data,
                                                char *jname,
                                                char *term_name);
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
    /* TODO Hier weiter machen*/
    
    /* the evaluated problem itself */
    socp->prb = (struct hhmpc_socp_prb*)malloc(sizeof(struct hhmpc_socp_prb));
    if (NULL == socp->prb) {return NULL;}
        socp->prb->q = socp->constant[HHMPC_Q_KL];
        socp->prb->r = socp->constant[HHMPC_R_KL];
    
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
    hhmpc_get_json_term(socp->constant[HHMPC_Q_KL], data, "constant", "q");
    hhmpc_get_json_term(socp->constant[HHMPC_R_KL], data, "constant", "r");
    
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
