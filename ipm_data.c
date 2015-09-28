#include "include/static_data.h"


uint32_t in_iter = 1;
uint32_t warm_start = 1;
real_t reg = 0.01;
// struct hhmpc_ipm_conf conf = {in_iter, warm_start, reg};
struct hhmpc_ipm_conf conf = {1, 1, 0.01};
struct hhmpc_ipm_P_hat P;
real_t kappa = 90.;
real_t roh = .1;

real_t z_opt[HHMPC_OS];
real_t delta_z[HHMPC_OS];
real_t v_opt[HHMPC_DS];
real_t delta_v[HHMPC_DS];

real_t h_hat[HHMPC_NB_IEQ];
real_t P_hat[HHMPC_NB_IEQ*HHMPC_OS];
real_t P_hat_T[HHMPC_OS*HHMPC_NB_IEQ];
real_t P2_hat[HHMPC_NB_IEQ*HHMPC_OS];
real_t P2_hat_T[HHMPC_OS*HHMPC_NB_IEQ];

real_t d[HHMPC_NB_IEQ];
real_t dsoft[HHMPC_NB_SOFT];
real_t diag_d_sq[HHMPC_NB_IEQ*HHMPC_NB_IEQ];
real_t diag_d_soft[HHMPC_NB_SOFT*HHMPC_NB_SOFT];
real_t Phi[HHMPC_OS*HHMPC_OS];
real_t Phi_soft[HHMPC_OS*HHMPC_OS];

real_t A_T[HHMPC_SV*HHMPC_SV];
real_t B_T[HHMPC_CV*HHMPC_SV];
real_t A_B_T[HHMPC_OV*HHMPC_SV];
real_t A_B[HHMPC_SV*HHMPC_OV];
real_t C_T[HHMPC_OS*HHMPC_DS];

real_t r_d[HHMPC_OS];
real_t r_d_soft[HHMPC_OS];
real_t r_p[HHMPC_DS];

real_t st_size[1];

real_t tmp1_optvar_seqlen[HHMPC_OS];
real_t tmp2_optvar_seqlen[HHMPC_OS];
real_t tmp2_dual_seqlen[HHMPC_DS];
real_t tmp3_mtx_optvar_nb_of_ueq[HHMPC_OS*HHMPC_NB_IEQ];
real_t tmp3_mtx_optvar_nb_of_soft[HHMPC_OS*HHMPC_NB_SOFT];
real_t tmp3_state_veclen[HHMPC_SV];
real_t tmp4_mtx_optvar_optvar[HHMPC_OS*HHMPC_OS];
real_t tmp4_nb_of_constr[HHMPC_NB_IEQ];
real_t tmp5_nb_of_constr[HHMPC_NB_IEQ];
real_t tmp6_optvar_seqlen[HHMPC_OS];
real_t tmp7_dual_seqlen[HHMPC_DS];
real_t tmp8_L_Y[(2*HHMPC_HORIZON-1)*HHMPC_SV*HHMPC_SV];
real_t tmp9_L_Y_T[(2*HHMPC_HORIZON-1)*HHMPC_SV*HHMPC_SV];
real_t tmp8_L_Phi[HHMPC_HORIZON*HHMPC_OV*HHMPC_OV];
real_t tmp9_L_Phi_T[HHMPC_HORIZON*HHMPC_OV*HHMPC_OV];
real_t tmp_phibl1[HHMPC_OV*HHMPC_OV];
real_t tmp_phibl2[HHMPC_OV*HHMPC_OV];
real_t tmp_phibl3[HHMPC_OV*HHMPC_OV];
real_t tmp10[HHMPC_OV*HHMPC_OV];
real_t tmpYbl[HHMPC_SV*HHMPC_SV];
real_t tmpQbl[HHMPC_SV*HHMPC_SV];
real_t tmp1_res_os[HHMPC_OS];
real_t tmp2_res_os[HHMPC_OS];
real_t tmp3_res_ds[HHMPC_DS];

real_t eye_optvar_veclen[HHMPC_OV*HHMPC_OV];
real_t eye_state_veclen[HHMPC_SV*HHMPC_SV];


struct hhmpc_ipm_socc *psocc[5];
struct hhmpc_ipm_socc socc[5];

real_t soccA_T1[30*30];
real_t soccA_T2[30*30];
real_t soccA_T3[30*30];
real_t soccA_T4[30*30];
real_t soccA_T5[30*30];
real_t *soccA_T[HHMPC_NB_SOCC];
real_t soccAAmcc1[30*30];
real_t soccAAmcc2[30*30];
real_t soccAAmcc3[30*30];
real_t soccAAmcc4[30*30];
real_t soccAAmcc5[30*30];
real_t *soccAAmcc[HHMPC_NB_SOCC];

void form_ipm(struct hhmpc_ipm *ipm, struct hhmpc_socp_prb *prb){
    
    uint32_t i;
    real_t *tmp;
    ipm->conf = &conf;
    ipm->P_of_z = &P;
    ipm->kappa = &kappa;
    ipm->roh = &roh;
    
    ipm->b = prb->b->data;
    /*ipm->h = prb->h->data;*/
    ipm->hsoft = prb->hsoft->data;
    ipm->Fusoft = prb->Fusoft->data;
    ipm->Fxsoft = prb->Fxsoft->data;
    ipm->Ffsoft = prb->Ffsoft->data;
    ipm->g = prb->g->data;
    ipm->A = prb->A->data;
    ipm->B = prb->B->data;
    ipm->C = prb->C->data;
    ipm->H = prb->H->data;
    ipm->z_ini = prb->z_ini->data;
    ipm->v_ini = prb->v_ini->data;
    ipm->zref = prb->zref->data;
    
    ipm->horizon = prb->horizon;
    ipm->state_veclen = prb->state_veclen;
    ipm->optvar_veclen = prb->optvar_veclen;
    ipm->optvar_seqlen = prb->optvar_seqlen;
    
    ipm->dual_seqlen = ipm->state_veclen * ipm->horizon;
    ipm->control_veclen = ipm->optvar_veclen - ipm->state_veclen;
    ipm->sizeof_optvar_seqlen = sizeof(real_t) * ipm->optvar_seqlen;
    ipm->sizeof_dual_seqlen = sizeof(real_t) * ipm->dual_seqlen;
    
    ipm->z_opt = z_opt;
    ipm->delta_z = delta_z;
    ipm->v_opt = v_opt;
    ipm->delta_v = delta_v;
    /* u_k points on first entries in z_opt*/
    prb->u_k->data = ipm->z_opt;
    prb->u_k->rows = ipm->control_veclen;
    prb->u_k->cols = 1;
    
    ipm->P_of_z->h = prb->h->data;
    ipm->P_of_z->P = prb->P->data;
    ipm->P_of_z->nb_lin_constr = prb->P->rows;
    ipm->P_of_z->nb_socc = prb->nb_socc;
    ipm->P_of_z->nb_qc = prb->nb_qc;
    ipm->nb_of_ueq_constr = prb->P->rows+prb->nb_socc+prb->nb_qc;
    ipm->nb_of_soft_constr =
            (prb->Fusoft->rows)*ipm->horizon + prb->Ffsoft->rows;
    ipm->rowsFusoft = prb->Fusoft->rows;
    ipm->rowsFfsoft = prb->Ffsoft->rows;
    
    ipm->Psoft = prb->Psoft->data;
    ipm->Psoft_T = prb->Psoft_T->data;
    mpcinc_mtx_transpose(ipm->Psoft_T, ipm->Psoft, prb->Psoft->rows, prb->Psoft->cols);
    ipm->P_of_z->h_hat = h_hat;
    ipm->P_of_z->socc = psocc;
    
    soccA_T[0] = soccA_T1;
    soccA_T[1] = soccA_T2;
    soccA_T[2] = soccA_T3;
    soccA_T[3] = soccA_T4;
    soccA_T[4] = soccA_T5;
    soccAAmcc[0] = soccAAmcc1;
    soccAAmcc[1] = soccAAmcc2;
    soccAAmcc[2] = soccAAmcc3;
    soccAAmcc[3] = soccAAmcc4;
    soccAAmcc[4] = soccAAmcc5;

    for (i = 0; i < ipm->P_of_z->nb_socc; i++){
        ipm->P_of_z->socc[i] = &socc[i];
        ipm->P_of_z->socc[i]->rowsA = prb->socc[i]->A->rows;
        ipm->P_of_z->socc[i]->colsA = prb->socc[i]->A->cols;
        ipm->P_of_z->socc[i]->A = prb->socc[i]->A->data;
        ipm->P_of_z->socc[i]->b = prb->socc[i]->b->data;
        ipm->P_of_z->socc[i]->c = prb->socc[i]->c->data;

        ipm->P_of_z->socc[i]->d = prb->socc[i]->d->data;
        ipm->P_of_z->socc[i]->par = ipm->z_opt+prb->socc[i]->par_0;
        ipm->P_of_z->socc[i]->par_0 = prb->socc[i]->par_0;
        ipm->P_of_z->socc[i]->par_l = prb->socc[i]->par_l;
        ipm->P_of_z->socc[i]->A_T = soccA_T[i];
        mpcinc_mtx_transpose(ipm->P_of_z->socc[i]->A_T, ipm->P_of_z->socc[i]->A,
                             prb->socc[i]->A->rows, prb->socc[i]->A->cols);
        ipm->P_of_z->socc[i]->AAmcc = soccAAmcc[i];
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
        tmp = ipm->P_of_z->h_hat+prb->P->rows+prb->nb_qc+i;
        mpcinc_mtx_multiply_mtx_vec(tmp, ipm->P_of_z->socc[i]->b, ipm->P_of_z->socc[i]->b,
                                    1, ipm->P_of_z->socc[i]->rowsA);
        tmp[0] = ipm->P_of_z->socc[i]->d[0]*ipm->P_of_z->socc[i]->d[0] - tmp[0];
    }
    
    ipm->P_of_z->P_hat = P_hat;
    zeroes(ipm->P_of_z->P_hat, ipm->optvar_seqlen*ipm->nb_of_ueq_constr);
    ipm->P_of_z->P_hat_T = P_hat_T;
    zeroes(ipm->P_of_z->P_hat_T, ipm->nb_of_ueq_constr*ipm->optvar_seqlen);
    ipm->P_of_z->P2_hat = P2_hat;
    zeroes(ipm->P_of_z->P2_hat, ipm->optvar_seqlen*ipm->nb_of_ueq_constr);
    ipm->P_of_z->P2_hat_T = P2_hat_T;
    zeroes(ipm->P_of_z->P2_hat_T, ipm->nb_of_ueq_constr*ipm->optvar_seqlen);
    
    ipm->h = ipm->P_of_z->h_hat;
    ipm->P = ipm->P_of_z->P_hat;
    ipm->P_T = ipm->P_of_z->P_hat_T;
    ipm->P2 = ipm->P_of_z->P2_hat;
    ipm->P2_T = ipm->P_of_z->P2_hat_T;
    
    ipm->j_in = &(ipm->conf->in_iter);
    ipm->reg = &(ipm->conf->reg);
    
    ipm->d = d;
    ipm->dsoft = dsoft;
    ipm->diag_d_sq = diag_d_sq;
    ipm->diag_d_soft = diag_d_soft;
    ipm->Phi = Phi;
    zeroes(ipm->Phi, ipm->optvar_seqlen*ipm->optvar_seqlen);
    ipm->Phi_soft = Phi_soft;
    zeroes(ipm->Phi_soft, ipm->optvar_seqlen*ipm->optvar_seqlen);
    ipm->A_T = A_T;
    mpcinc_mtx_transpose(ipm->A_T, ipm->A, prb->A->rows, prb->A->cols);
    ipm->B_T = B_T;
    mpcinc_mtx_transpose(ipm->B_T, ipm->B, prb->B->rows, prb->B->cols);
    ipm->A_B_T = A_B_T;
    for (i = 0; i < ipm->state_veclen*ipm->state_veclen; i++)
        ipm->A_B_T[i] = ipm->A_T[i];
    for (i = 0; i < ipm->state_veclen*ipm->control_veclen; i++)
        ipm->A_B_T[ipm->state_veclen*ipm->state_veclen+i] = ipm->B_T[i];
    ipm->A_B = A_B;
    mpcinc_mtx_transpose(ipm->A_B, ipm->A_B_T, ipm->optvar_veclen, ipm->state_veclen);
    ipm->C_T = C_T;
    mpcinc_mtx_transpose(ipm->C_T, ipm->C, prb->C->rows, prb->C->cols);
    
    ipm->r_d = r_d;
    ipm->r_d_soft = r_d_soft;
    ipm->r_p = r_p;
    ipm->st_size = st_size;
    
    ipm->tmp1_optvar_seqlen = tmp1_optvar_seqlen;
    ipm->tmp2_optvar_seqlen = tmp2_optvar_seqlen;
    ipm->tmp2_dual_seqlen = tmp2_dual_seqlen;
    ipm->tmp3_mtx_optvar_nb_of_ueq = tmp3_mtx_optvar_nb_of_ueq;
    ipm->tmp3_mtx_optvar_nb_of_soft = tmp3_mtx_optvar_nb_of_soft;
    ipm->tmp3_state_veclen = tmp3_state_veclen;
    ipm->tmp4_mtx_optvar_optvar = tmp4_mtx_optvar_optvar;
    ipm->tmp4_nb_of_constr = tmp4_nb_of_constr;
    ipm->tmp5_nb_of_constr = tmp5_nb_of_constr;
    ipm->tmp6_optvar_seqlen = tmp6_optvar_seqlen;
    ipm->tmp7_dual_seqlen = tmp7_dual_seqlen;
    ipm->tmp8_L_Y = tmp8_L_Y;
    ipm->tmp9_L_Y_T = tmp9_L_Y_T;
    ipm-> tmp8_L_Phi= tmp8_L_Phi;
    ipm->tmp9_L_Phi_T = tmp9_L_Phi_T;
    ipm->tmp_phibl1 = tmp_phibl1;
    ipm->tmp_phibl2 = tmp_phibl2;
    ipm->tmp_phibl3 = tmp_phibl3;
    ipm->tmp10 = tmp10;
    ipm->tmpYbl = tmpYbl;
    ipm->tmpQbl = tmpQbl;
    ipm->tmp1_res_os = tmp1_res_os;
    ipm->tmp2_res_os = tmp2_res_os;
    ipm->tmp3_res_ds = tmp3_res_ds;

    ipm->eye_optvar_veclen = eye_optvar_veclen;
    eye(ipm->eye_optvar_veclen, ipm->optvar_veclen);
    ipm->eye_state_veclen = eye_state_veclen;
    eye(ipm->eye_state_veclen, ipm->state_veclen);
    
}