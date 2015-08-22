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
    real_t *q;
    real_t *r;
};


#endif /* HHMPCIPM_H */