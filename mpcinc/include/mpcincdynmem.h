#ifndef MPCINC_DYNMEM_H
#define MPCINC_DYNMEM_H

#include "cjson.h"

typedef enum mpcinc_dynmem_error {
    MPCINC_DYNMEM_OK = 0,
    MPCINC_DYNMEM_FAIL = 1,
} mpcinc_dynmem_error_t;

cJSON *mpcinc_dynmem_get_data(char *fname);

#endif
