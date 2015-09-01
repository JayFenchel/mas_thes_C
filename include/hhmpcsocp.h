#ifndef HHMPCSOCP_H
#define HHMPCSOCP_H

#define HIER    printf("HIER\n");

#include "mc04types.h"  /* typedefs */
#include "arithmetic.h"

/*Festlegen, was Parameter, Konstante oder parametrisch ist*/
enum {
        HHMPC_XK, /* State of the system. */
    
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
        HHMPC_H,
        HHMPC_C,
        
    HHMPC_CONST_NUM
};

enum {
    
    HHMPC_B_KL,
    HHMPC_H_KL,
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
    struct hhmpc_term *xr;
    struct hhmpc_term *x_k;
    struct hhmpc_term *ur;
    struct hhmpc_term *u_lb;
    struct hhmpc_term *u_ub;
    /*eigene*/
    struct hhmpc_term *g;
    struct hhmpc_term *q;
    struct hhmpc_term *r;
    struct hhmpc_term *S;
    struct hhmpc_term *S_T;
    
    struct hhmpc_term *A;
    struct hhmpc_term *B;
    struct hhmpc_term *C;
    struct hhmpc_term *P;
    uint32_t nb_qc;
    struct hhmpc_qc **qc;
    uint32_t nb_socc;
    struct hhmpc_socc **socc;
    struct hhmpc_term *H;
    struct hhmpc_term *b;
    
    struct hhmpc_term *h;
    /* To calculate length of vectors such as g */
    uint32_t horizon; /* Prediction horizon. */
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


#endif /* HHMPCSOCP_H */