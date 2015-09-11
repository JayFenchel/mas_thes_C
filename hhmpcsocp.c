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

void sim_next_xk(const struct hhmpc_socp *socp)
{   
    struct hhmpc_term *tmp = socp->prb->tmp_state_veclen;
    hhmpc_copy_data(tmp, socp->prb->x_k);
    mpcinc_mtx_multiply_mtx_vec(socp->prb->x_k->data,
                                socp->prb->A->data, tmp->data,
                                socp->prb->A->rows, socp->prb->A->cols);
    mpcinc_mtx_multiply_mtx_vec(tmp->data,
                                socp->prb->B->data, socp->prb->u_k->data,
                                socp->prb->B->rows, socp->prb->B->cols);
    mpcinc_mtx_add_direct(socp->prb->x_k->data, tmp->data, socp->prb->x_k->rows, 1);
}

void hhmpc_copy_data(struct hhmpc_term *dest, struct hhmpc_term *src)  {
    int j;
    for (j=0; j<(dest->cols*dest->rows); j++) {
        dest->data[j] = src->data[j];
    }
    return;
}

