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
    real_t *q;
    real_t *r;
    
    uint32_t *j_in;
    
    real_t *kappa;  /* Barrier parameter. */
    uint32_t horizon; /* Prediction horizon. */
    uint32_t optvar_veclen;  /* The length of each vector in the optimation variable sequence. */
    uint32_t optvar_seqlen;  /* The full length of optimization variable sequence. */
    uint32_t sizeof_optvar_seqlen;  /* Number of bytes in the optimization variable sequence. */
    
    real_t *tmp1_optvar_seqlen;  /* Temporary variable of length optvar_seqlen. */
};


#endif /* HHMPCIPM_H */