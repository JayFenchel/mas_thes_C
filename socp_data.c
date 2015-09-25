#include "include/static_data.h"
#include "include/hhmpcusefull.h"
// #include <stdio.h>

// #include "hhmpc_testdata.h"
#ifdef HHMPC_QPTEST
#include "hhmpc_qpdata01.h"
#endif
#ifdef HHMPC_QPSMALLTEST
#include "hhmpc_qpdatasmall02.h"
#endif
#ifdef HHMPC_SOCPTEST
#include "hhmpc_socpdata01.h"
#endif
#ifdef HHMPC_SOCPHARDTEST
#include "hhmpc_socpharddata01.h"
#endif
#ifdef HHMPC_SOCPSOFTTEST
#include "hhmpc_socpwithsoftdata01.h"
#endif



void form_socp_H(struct hhmpc_socp *socp)
{
    socp->constant[HHMPC_H] = &H_term;
}
void form_socp(struct hhmpc_socp *socp)
{

socp->constant[HHMPC_Q_KL] = &q_term;
socp->constant[HHMPC_Q] = &Q_term;
socp->constant[HHMPC_R_KL] = &r_term;
socp->constant[HHMPC_R] = &R_term;
socp->constant[HHMPC_S] = &S_term;
socp->constant[HHMPC_S_T] = &S_T_term;  /* 2*S_T */
socp->constant[HHMPC_A] = &A_term;
socp->constant[HHMPC_B] = &B_term;
socp->constant[HHMPC_PSOFT] = &Psoft_term;
socp->constant[HHMPC_P] = &P_term;
socp->constant[HHMPC_FXSOFT] = &Fxsoft_term;
socp->constant[HHMPC_FUSOFT] = &Fusoft_term;
socp->constant[HHMPC_FSOFT] = &fsoft_term;
socp->constant[HHMPC_FFSOFT] = &Ffsoft_term;
socp->constant[HHMPC_FFSOFT_KL] = &ffsoft_term;
socp->constant[HHMPC_H] = &H_term;

socp->constant[HHMPC_C] = &C_term;

socp->par[HHMPC_XK] = &xk_term;
socp->par[HHMPC_XR] = &xr_term;
socp->par[HHMPC_UR] = &ur_term;
socp->par[HHMPC_ZINI] = &zini_term;
socp->par[HHMPC_VINI] = &vini_term;

socp->pmetric[HHMPC_ZR] = &zref;
socp->pmetric[HHMPC_ZR]->fac_num = &zref_fac_num;
socp->pmetric[HHMPC_ZR]->val = &zref_val_term;
socp->pmetric[HHMPC_ZR]->aux = &zref_aux_term;
socp->pmetric[HHMPC_ZR]->fac0 = &zref_fac0_term;
socp->pmetric[HHMPC_ZR]->par = pzref_par_term;
socp->pmetric[HHMPC_ZR]->fac = pzref_fac_term;
socp->pmetric[HHMPC_ZR]->fac[0] = &zref_fac1_term;
socp->pmetric[HHMPC_ZR]->par[0] = &ur_term;
socp->pmetric[HHMPC_ZR]->fac[1] = &zref_fac2_term;
socp->pmetric[HHMPC_ZR]->par[1] = &xr_term;

socp->pmetric[HHMPC_B_KL] = &b;
socp->pmetric[HHMPC_B_KL]->fac_num = &b_fac_num;
socp->pmetric[HHMPC_B_KL]->val = &b_val_term;
socp->pmetric[HHMPC_B_KL]->aux = &b_aux_term;
socp->pmetric[HHMPC_B_KL]->fac0 = &b_fac0_term;
socp->pmetric[HHMPC_B_KL]->par = &pb_par_term;
socp->pmetric[HHMPC_B_KL]->fac = &pb_fac_term;
socp->pmetric[HHMPC_B_KL]->fac[0] = &b_fac1_term;
socp->pmetric[HHMPC_B_KL]->par[0] = &xk_term;

socp->pmetric[HHMPC_H_KL] = &h;
socp->pmetric[HHMPC_H_KL]->fac_num = &h_fac_num;
socp->pmetric[HHMPC_H_KL]->val = &h_val_term;
socp->pmetric[HHMPC_H_KL]->aux = &h_aux_term;
socp->pmetric[HHMPC_H_KL]->fac0 = &h_fac0_term;
socp->pmetric[HHMPC_H_KL]->par = &ph_par_term;
socp->pmetric[HHMPC_H_KL]->fac = &ph_fac_term;
socp->pmetric[HHMPC_H_KL]->fac[0] = &h_fac1_term;
socp->pmetric[HHMPC_H_KL]->par[0] = &xk_term;

socp->pmetric[HHMPC_HSOFT_KL] = &hsoft;
socp->pmetric[HHMPC_HSOFT_KL]->fac_num = &hsoft_fac_num;
socp->pmetric[HHMPC_HSOFT_KL]->val = &hsoft_val_term;
socp->pmetric[HHMPC_HSOFT_KL]->aux = &hsoft_aux_term;
socp->pmetric[HHMPC_HSOFT_KL]->fac0 = &hsoft_fac0_term;
socp->pmetric[HHMPC_HSOFT_KL]->par = &phsoft_par_term;
socp->pmetric[HHMPC_HSOFT_KL]->fac = &phsoft_fac_term;
socp->pmetric[HHMPC_HSOFT_KL]->fac[0] = &hsoft_fac1_term;
socp->pmetric[HHMPC_HSOFT_KL]->par[0] = &xk_term;

socp->pmetric[HHMPC_G_KL] = &g;
socp->pmetric[HHMPC_G_KL]->fac_num = &g_fac_num;
socp->pmetric[HHMPC_G_KL]->val = &g_val_term;
socp->pmetric[HHMPC_G_KL]->aux = &g_aux_term;
socp->pmetric[HHMPC_G_KL]->fac0 = &g_fac0_term;
socp->pmetric[HHMPC_G_KL]->par = &pg_par_term;
socp->pmetric[HHMPC_G_KL]->fac = &pg_fac_term;
socp->pmetric[HHMPC_G_KL]->fac[0] = socp->constant[HHMPC_S_T];
socp->pmetric[HHMPC_G_KL]->par[0] = &xk_term;
socp->prb = &prb;

socp->prb->z_ini = socp->par[HHMPC_ZINI];
socp->prb->v_ini = socp->par[HHMPC_VINI];
socp->prb->x_k = socp->par[HHMPC_XK];
socp->prb->tmp_state_veclen = &tmp_state_veclen_term;

socp->prb->zref = socp->pmetric[HHMPC_ZR]->val;
socp->prb->b = socp->pmetric[HHMPC_B_KL]->val;
socp->prb->h = socp->pmetric[HHMPC_H_KL]->val;
socp->prb->hsoft = socp->pmetric[HHMPC_HSOFT_KL]->val;
socp->prb->Fusoft = socp->constant[HHMPC_FUSOFT];
socp->prb->Fxsoft = socp->constant[HHMPC_FXSOFT];
socp->prb->Ffsoft = socp->constant[HHMPC_FFSOFT];
socp->prb->Psoft = socp->constant[HHMPC_PSOFT];        
socp->prb->Psoft_T = &Psoft_T_term;        
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

socp->prb->socc = psocc_term;
socp->prb->socc[0] = &socc1_term; 
socp->prb->socc[1] = &socc2_term; 
socp->prb->socc[2] = &socc3_term; 
socp->prb->socc[3] = &socc4_term; 
socp->prb->socc[4] = &socc5_term; 

socp->prb->u_k = &u_k_term;

socp->prb->nb_qc = nb_qc;
socp->prb->nb_socc = nb_socc;
socp->prb->horizon = horizon;
socp->prb->state_veclen = state_veclen;
socp->prb->optvar_veclen = optvar_veclen;
socp->prb->optvar_seqlen = optvar_seqlen;
sizeof_optvar_seqlen = sizeof(real_t) * optvar_seqlen; //eigenttlich kein sizeof verwenden
socp->prb->sizeof_optvar_seqlen = sizeof_optvar_seqlen;

}
