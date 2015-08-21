#ifndef HHMPCSOCP_H
#define HHMPCSOCP_H

#include "mc04types.h"  /* typedefs */
#include "arithmetic.h"


struct hhmpc_term {
    uint32_t rows;
    uint32_t cols;
    real_t *data;
};

struct hhmpc_socp_prb {
    struct hhmpc_term *xr;
    struct hhmpc_term *x_k;
    struct hhmpc_term *ur;
    struct hhmpc_term *g;
    struct hhmpc_term *H;
    struct hhmpc_term *u_lb;
    struct hhmpc_term *u_ub;
};

struct hhmpc_socp {
    /*struct hhmpc_term *par[MPCINC_PAR_NUM];
    struct hhmpc_term *constant[MPCINC_CONSTANT_NUM];
    struct mpcinc_pmetric *pmetric[MPCINC_PMETRIC_NUM];*/
    struct hhmpc_socp_prb *prb;
};


#endif /* HHMPCSOCP_H */