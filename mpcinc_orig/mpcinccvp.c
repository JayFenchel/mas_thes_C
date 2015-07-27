#include "mpcincmtxops.h"
#include "mpcinccvp.h"


static void mpcinc_copy_data(struct mpcinc_term *dest, struct mpcinc_term *src);

void mpcinc_cvp_form_problem(struct mpcinc_cvp *cvp)  {
  int i, j;
  struct mpcinc_pmetric *pm;

    for (i=0; i<MPCINC_PMETRIC_NUM; i++) {
        pm = cvp->pmetric[i];
        mpcinc_copy_data(pm->val, pm->fac0);
        for (j=0; j<pm->fac_num[0]; j++) {
            mpcinc_mtx_mul_add(pm->val->data, pm->aux->data,
            pm->fac[j]->data, pm->par[j]->data,
            pm->fac[j]->rows, pm->fac[j]->cols);
        }
    }
    return;
}

void mpcinc_copy_data(struct mpcinc_term *dest, struct mpcinc_term *src)  {
    int j;
    for (j=0; j<(dest->cols*dest->rows); j++) {
        dest->data[j] = src->data[j];
    }
    return;
}

