#ifndef HHMPCSOCP_H
#define HHMPCSOCP_H

#define HIER    printf("HIER\n");

#include "mc04types.h"  /* typedefs */
#include "arithmetic.h"

/*Festlegen, was Parameter, Konstante oder parametrisch ist*/
enum {
        HHMPC_XK, /* State of the system. */
        HHMPC_XR, /* Reference for states */
        HHMPC_UR, /* Reference for controls */
        HHMPC_ZINI,  /* Initial guess for opt var z */
        HHMPC_VINI,  /* Initial guess for dual var v */
    
    HHMPC_PAR_NUM
};

enum {
        HHMPC_Q_KL,
        HHMPC_Q,
        HHMPC_R_KL,
        HHMPC_R,
        HHMPC_S,
        HHMPC_S_T,
        HHMPC_A,
        HHMPC_B,
        
        HHMPC_P,
        /* Soft Constraints: */
        HHMPC_PSOFT,
        HHMPC_FXSOFT,  /* mtx[Fx, Fu]*[x, u] <= f */
        HHMPC_FUSOFT,
        HHMPC_FSOFT,
        HHMPC_FFSOFT,  /* Ff*x(t+T) <= ff */
        HHMPC_FFSOFT_KL,
        HHMPC_H,
        HHMPC_C,
        
    HHMPC_CONST_NUM
};

enum {
    HHMPC_ZR,
    HHMPC_B_KL,
    HHMPC_H_KL,
    HHMPC_HSOFT_KL,
    HHMPC_G_KL,
    
    HHMPC_PMETRIC_NUM
};


struct hhmpc_term {
    uint32_t rows;
    uint32_t cols;
    real_t *data;
};

struct hhmpc_pmetric {
    uint32_t *fac_num;
    struct hhmpc_term *val;
    struct hhmpc_term *aux;
    struct hhmpc_term *fac0;
    struct hhmpc_term **fac;
    struct hhmpc_term **par;
};

/**/
struct hhmpc_qc {
    struct hhmpc_term *Gamma;
    struct hhmpc_term *beta;
    struct hhmpc_term *alpha;
    uint32_t par_0;
    uint32_t par_l;
};

/**/
struct hhmpc_socc {
    struct hhmpc_term *A;
    struct hhmpc_term *b;
    struct hhmpc_term *c;
    struct hhmpc_term *d;
    uint32_t par_0;
    uint32_t par_l;
};

struct hhmpc_socp_prb {
    struct hhmpc_term *x_k;
    struct hhmpc_term *tmp_state_veclen;
    struct hhmpc_term *u_k;
    struct hhmpc_term *zref;
    struct hhmpc_term *z_ini;
    struct hhmpc_term *v_ini;
    
    struct hhmpc_term *g;
    struct hhmpc_term *q;
    struct hhmpc_term *r;
    struct hhmpc_term *S;
    struct hhmpc_term *S_T;
    
    struct hhmpc_term *A;
    struct hhmpc_term *B;
    struct hhmpc_term *C;
    struct hhmpc_term *P;
    struct hhmpc_term *Psoft;
    struct hhmpc_term *Psoft_T;
    struct hhmpc_term *Fusoft;
    struct hhmpc_term *Fxsoft;
    struct hhmpc_term *Ffsoft;
    struct hhmpc_term *fsoft;
    struct hhmpc_term *ffsoft;
    uint32_t nb_qc;
    struct hhmpc_qc **qc;
    uint32_t nb_socc;
    struct hhmpc_socc **socc;
    struct hhmpc_term *H;
    struct hhmpc_term *b;
    
    struct hhmpc_term *h;
    struct hhmpc_term *hsoft;
    /* To calculate length of vectors such as g */
    uint32_t horizon; /* Prediction horizon. */
    uint32_t state_veclen;
    uint32_t optvar_veclen;  /* The length of each vector in the optimation variable sequence. */
    uint32_t optvar_seqlen;  /* The full length of optimization variable sequence. */
    uint32_t sizeof_optvar_seqlen;  /* Number of bytes in the optimization variable sequence. */
};

struct hhmpc_socp {
    struct hhmpc_term *par[HHMPC_PAR_NUM];
    struct hhmpc_term *constant[HHMPC_CONST_NUM];
    struct hhmpc_pmetric *pmetric[HHMPC_PMETRIC_NUM];
    struct hhmpc_socp_prb *prb;
};

extern void hhmpc_socp_form_problem(struct hhmpc_socp *socp);

extern void sim_next_xk(const struct hhmpc_socp *socp);


#endif /* HHMPCSOCP_H */