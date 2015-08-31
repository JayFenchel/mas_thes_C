#ifndef HHMPCIPM_H
#define HHMPCIPM_H


#include "mc04types.h"
#include "arithmetic.h"

/* Configuration parameters of the HHMPC algorithm. */
struct hhmpc_ipm_conf {
    uint32_t in_iter;
};

struct hhmpc_ipm {
    struct hhmpc_ipm_conf *conf;  /* Algorithm configuration data. */
    real_t *z_opt;  /* Solution to the optimal control problem. */
    real_t *z_ini;  /* Initial guess for the optimal control sequence. */
    real_t *delta_z;
    real_t *v_ini;
    real_t *v_opt;
    real_t *delta_v;
    real_t *q;
    real_t *r;
    real_t *g;
    real_t *A;
    real_t *A_T;
    real_t *B;
    real_t *B_T;
    real_t *C;
    real_t *C_T;
    real_t *H;
    real_t *P;
    real_t *P_T;
    
    real_t *b;
    real_t *h;
    
    real_t *d;
    real_t *diag_d_sq;
    real_t *Phi;
    real_t *r_p;
    real_t *r_d;
    
    real_t *st_size;
    uint32_t *j_in;
    
    real_t *kappa;  /* Barrier parameter. */
    uint32_t horizon; /* Prediction horizon. */
    uint32_t optvar_veclen;  /* The length of each vector in the optimation variable sequence. */
    uint32_t optvar_seqlen;  /* The full length of optimization variable sequence. */
    uint32_t state_veclen;  /* Dimension of state variable x_k. */
    uint32_t dual_seqlen;  /* Full length of dual variable v associated with the eq constr. */
    uint32_t nb_of_ueq_constr;
    uint32_t sizeof_dual_seqlen;
    uint32_t sizeof_optvar_seqlen;  /* Number of bytes in the optimization variable sequence. */
    
    real_t *tmp1_optvar_seqlen;  /* Temporary variable of length optvar_seqlen. */
    real_t *tmp2_optvar_seqlen;  /* Temporary variable of length optvar_seqlen. */
    real_t *tmp2_dual_seqlen;  /* Temporary variable of length optvar_seqlen. */
    real_t *tmp3_state_veclen;
    real_t *tmp3_mtx_optvar_nb_of_ueq;
    real_t *tmp4_nb_of_constr;
    real_t *tmp5_nb_of_constr;
    real_t *tmp6_optvar_seqlen;  /* Temporary variable of length optvar_seqlen. */
    real_t *tmp7_dual_seqlen;  /* Temporary variable of length optvar_seqlen. */
    real_t *tmp8_L_Y;
    real_t *tmp9_L_Y_T;
    real_t *eye_optvar_veclen;
    real_t *eye_state_veclen;
};

/* External function declarations */

extern void hhmpc_ipm_solve_problem(const struct hhmpc_ipm *ipm);

extern uint32_t hhmpc_ipm_check_valid(const struct hhmpc_ipm *ipm, const real_t *z_check);


#endif /* HHMPCIPM_H */