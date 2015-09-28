#ifndef MPCFORMQP_H
#define MPCFORMQP_H

#include "mc04types.h"  /* typedefs */
#include "arithmetic.h"

/* Declarations that depend only on the structure of the problem
 * These are used to generate code independent of the data size
 */
enum {
        MPC_X_K,

    MPC_PAR_NUM
};

enum {
        MPC_U_LB,
        MPC_WVT4,
        MPC_WM3,
        MPC_H,
        MPC_WM1,
        MPC_WM0,
        MPC_WM4,
        MPC_WVT3,
        MPC_WVT2,
        MPC_WVT1,
        MPC_WVT0,
        MPC_WM2,
        MPC_U_UB,

    MPC_CONSTANT_NUM
};

enum {
        MPC_WN2,
        MPC_G,
        MPC_WN3,
        MPC_WS3,
        MPC_WS1,
        MPC_WS0,
        MPC_WN0,
        MPC_WN1,
        MPC_WS4,
        MPC_WN4,
        MPC_WS2,

    MPC_PMETRIC_NUM
};

struct mpc_term {
    uint32_t rows;
    uint32_t cols;
    real_t *data;
};

struct mpc_pmetric {
    uint32_t *fac_num;
    struct mpc_term *val;
    struct mpc_term *aux;
    struct mpc_term *fac0;
    struct mpc_term **fac;
    struct mpc_term **par;
};

struct mpc_cvp_prb {
    struct mpc_term *wn2;
    struct mpc_term *g;
    struct mpc_term *wn3;
    struct mpc_term *ws3;
    struct mpc_term *ws1;
    struct mpc_term *ws0;
    struct mpc_term *wn0;
    struct mpc_term *wn1;
    struct mpc_term *ws4;
    struct mpc_term *wn4;
    struct mpc_term *ws2;
    struct mpc_term *u_lb;
    struct mpc_term *wvT4;
    struct mpc_term *Wm3;
    struct mpc_term *H;
    struct mpc_term *Wm1;
    struct mpc_term *Wm0;
    struct mpc_term *Wm4;
    struct mpc_term *wvT3;
    struct mpc_term *wvT2;
    struct mpc_term *wvT1;
    struct mpc_term *wvT0;
    struct mpc_term *Wm2;
    struct mpc_term *u_ub;

#if 0  /* left as reference */
    struct mpc_term *H;  /**< The Hessian matrix. */
    struct mpc_term *g;  /**< The gradient vector. */
    struct mpc_term *u_lb;  /**< The lower bound for the box constraints. */
    struct mpc_term *u_ub;  /**< The upper bound for the box constraints. */
    struct mpc_term *V;  /**< The mixed constraints matrix. */
    struct mpc_term *v_lb;  /**< The lower bound for the mixed constraints. */
    struct mpc_term *v_ub;  /**< The upper bound for the mixed constraints. */
#endif
};  /**< MPC quadratic program form for a given system state x. 
 * The quadratic program to solve has the form:
 * minimize 0.5 * u^T * HoL * u + u^T * goL
 * subject to u_lb <= u <= u_ub
 *        v_lb <= V * x <= v_ub
 *
 * the transpose of a matrix is denote by ^T.
 */

struct mpc_cvp {
    struct mpc_term *par[MPC_PAR_NUM];
    struct mpc_term *constant[MPC_CONSTANT_NUM];
    struct mpc_pmetric *pmetric[MPC_PMETRIC_NUM];
    struct mpc_cvp_prb *prb;
};

struct mpc_cvp_parameters {
    real_t *x_k;

};

extern void mpc_cvp_copy_parameters(struct mpc_cvp *cvp, struct mpc_cvp_parameters *p);
extern void mpc_cvp_form_problem(struct mpc_cvp *cvp);

#endif
