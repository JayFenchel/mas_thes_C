#include "include/mpcinccvp.h"
    real_t mpcinc_xr_data[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    struct mpcinc_term mpcinc_xr = {6, 1, mpcinc_xr_data};
    real_t mpcinc_x_k_data[] = {0.0, 0.0};
    struct mpcinc_term mpcinc_x_k = {2, 1, mpcinc_x_k_data};
    real_t mpcinc_ur_data[] = {0.0, 0.0};
    struct mpcinc_term mpcinc_ur = {2, 1, mpcinc_ur_data};
    real_t mpcinc_g_val_data[] = {0.0, 0.0};
    struct mpcinc_term mpcinc_g_val = {2, 1, mpcinc_g_val_data};
    real_t mpcinc_g_fac0_data[] = {0.0, 0.0};
    struct mpcinc_term mpcinc_g_fac0 = {2, 1, mpcinc_g_fac0_data};
    real_t mpcinc_g_aux_data[] = {0.0, 0.0};
    struct mpcinc_term mpcinc_g_aux = {2, 1, mpcinc_g_aux_data};
    struct mpcinc_term *mpcinc_g_fac[3];
    struct mpcinc_term *mpcinc_g_par[3];
    struct mpcinc_pmetric mpcinc_g;
    real_t mpcinc_g_xr_data[] = {0.0, 0.0, -9.6473220735445597e-05, -0.018978390346516102, -0.05207122612571069, -0.069520660328755937, 0.0, 0.0, 0.0, 0.0, -0.057128117481019196, -0.076202981048270751};
    struct mpcinc_term mpcinc_g_xr = {2, 6, mpcinc_g_xr_data};
    real_t mpcinc_g_x_k_data[] = {0.052167699346446134, 0.075035869563750038, 0.057128117481019196, 0.06342528114382355};
    struct mpcinc_term mpcinc_g_x_k = {2, 2, mpcinc_g_x_k_data};
    real_t mpcinc_g_ur_data[] = {-0.99715612769352779, 0.0, 0.0, -0.99715612769352779};
    struct mpcinc_term mpcinc_g_ur = {2, 2, mpcinc_g_ur_data};
    uint32_t mpcinc_g_fac_num = 3;
    real_t mpcinc_H_data[] = {0.9987290508536623, 0.001328190913114191, 0.0013281909131141912, 0.9986119892312272};
    struct mpcinc_term mpcinc_H = {2, 2, mpcinc_H_data};
    real_t mpcinc_u_lb_data[] = {-100.0, -100.0};
    struct mpcinc_term mpcinc_u_lb = {2, 1, mpcinc_u_lb_data};
    real_t mpcinc_u_ub_data[] = {100.0, 100.0};
    struct mpcinc_term mpcinc_u_ub = {2, 1, mpcinc_u_ub_data};
    struct mpcinc_cvp_prb mpcinc_prb;

void mpcinc_initialize_problem_structure(struct mpcinc_cvp *cvp) {
    cvp->par[MPCINC_XR] = &mpcinc_xr;
    cvp->par[MPCINC_X_K] = &mpcinc_x_k;
    cvp->par[MPCINC_UR] = &mpcinc_ur;
    cvp->pmetric[MPCINC_G] = &mpcinc_g;
    cvp->pmetric[MPCINC_G]->fac_num = &mpcinc_g_fac_num;
    cvp->pmetric[MPCINC_G]->val = &mpcinc_g_val;
    cvp->pmetric[MPCINC_G]->aux = &mpcinc_g_aux;
    cvp->pmetric[MPCINC_G]->fac0 = &mpcinc_g_fac0;
    cvp->pmetric[MPCINC_G]->par = mpcinc_g_par;
    cvp->pmetric[MPCINC_G]->fac = mpcinc_g_fac;
    cvp->pmetric[MPCINC_G]->fac[0] = &mpcinc_g_xr;
    cvp->pmetric[MPCINC_G]->par[0] = &mpcinc_xr;
    cvp->pmetric[MPCINC_G]->fac[1] = &mpcinc_g_x_k;
    cvp->pmetric[MPCINC_G]->par[1] = &mpcinc_x_k;
    cvp->pmetric[MPCINC_G]->fac[2] = &mpcinc_g_ur;
    cvp->pmetric[MPCINC_G]->par[2] = &mpcinc_ur;
    cvp->constant[MPCINC_H] = &mpcinc_H;
    cvp->constant[MPCINC_U_LB] = &mpcinc_u_lb;
    cvp->constant[MPCINC_U_UB] = &mpcinc_u_ub;

    cvp->prb = &mpcinc_prb;
    cvp->prb->xr = &mpcinc_xr;
    cvp->prb->x_k = &mpcinc_x_k;
    cvp->prb->ur = &mpcinc_ur;
    cvp->prb->g = &mpcinc_g_val;
    cvp->prb->H = &mpcinc_H;
    cvp->prb->u_lb = &mpcinc_u_lb;
    cvp->prb->u_ub = &mpcinc_u_ub;

    return;
}
