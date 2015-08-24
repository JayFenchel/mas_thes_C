#include "include/hhmpcsocp.h"
#include "include/mpcincmtxops.h"


static void hhmpc_copy_data(struct hhmpc_term *dest, struct hhmpc_term *src);

void hhmpc_socp_form_problem(struct hhmpc_socp *socp)
{
    /* Needed if parametric Matrices are used*/
    int i, j;
    struct hhmpc_pmetric *pm;
    
    for (i = 0;  i < HHMPC_PMETRIC_NUM; i++) {
        pm = socp->pmetric[i];
        hhmpc_copy_data(pm->val, pm->fac0);
        for (j = 0; j < pm->fac_num[0]; j++){
            mpcinc_mtx_mul_add(pm->val->data, pm->aux->data,
            pm->fac[j]->data, pm->par[j]->data,
            pm->fac[j]->rows, pm->fac[j]->cols);
        }
    }
    
}

void hhmpc_copy_data(struct hhmpc_term *dest, struct hhmpc_term *src)  {
    int j;
    for (j=0; j<(dest->cols*dest->rows); j++) {
        dest->data[j] = src->data[j];
    }
    return;
}

