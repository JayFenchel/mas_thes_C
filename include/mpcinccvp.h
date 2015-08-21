#ifndef MPCINCFORMQP_H
#define MPCINCFORMQP_H

#include "mc04types.h"  /* typedefs */
#include "arithmetic.h"

/* Declarations that depend only on the structure of the problem
 * These are used to generate code independent of the data size
 */
enum {
        MPCINC_XR,
        MPCINC_X_K,
        MPCINC_UR,

    MPCINC_PAR_NUM
};

enum {
        MPCINC_H,
        MPCINC_U_LB,
        MPCINC_U_UB,

    MPCINC_CONSTANT_NUM
};

enum {
        MPCINC_G,

    MPCINC_PMETRIC_NUM
};

struct mpcinc_term {
    uint32_t rows;
    uint32_t cols;
    real_t *data;
};

struct mpcinc_pmetric {
    uint32_t *fac_num;
    struct mpcinc_term *val;
    struct mpcinc_term *aux;
    struct mpcinc_term *fac0;
    struct mpcinc_term **fac;
    struct mpcinc_term **par;
};

struct mpcinc_cvp_prb {
    struct mpcinc_term *xr;
    struct mpcinc_term *x_k;
    struct mpcinc_term *ur;
    struct mpcinc_term *g;
    struct mpcinc_term *H;
    struct mpcinc_term *u_lb;
    struct mpcinc_term *u_ub;

#if 0  /* left as reference */
    struct mpcinc_term *H;  /**< The Hessian matrix. */
    struct mpcinc_term *g;  /**< The gradient vector. */
    struct mpcinc_term *u_lb;  /**< The lower bound for the box constraints. */
    struct mpcinc_term *u_ub;  /**< The upper bound for the box constraints. */
    struct mpcinc_term *V;  /**< The mixed constraints matrix. */
    struct mpcinc_term *v_lb;  /**< The lower bound for the mixed constraints. */
    struct mpcinc_term *v_ub;  /**< The upper bound for the mixed constraints. */
#endif
};  /**< MPCINC quadratic program form for a given system state x. 
 * The quadratic program to solve has the form:
 * minimize 0.5 * u^T * HoL * u + u^T * goL
 * subject to u_lb <= u <= u_ub
 *        v_lb <= V * x <= v_ub
 *
 * the transpose of a matrix is denote by ^T.
 */

struct mpcinc_cvp {
    struct mpcinc_term *par[MPCINC_PAR_NUM];
    struct mpcinc_term *constant[MPCINC_CONSTANT_NUM];
    struct mpcinc_pmetric *pmetric[MPCINC_PMETRIC_NUM];
    struct mpcinc_cvp_prb *prb;
};

extern void mpcinc_cvp_form_problem(struct mpcinc_cvp *cvp);

#endif
