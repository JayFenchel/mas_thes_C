#ifndef HHMPCDYNMEM_H
#define HHMPCDYNMEM_H


#include "cjson.h"

typedef enum hhmpc_dynmem_error {
    HHMPC_DYNMEM_OK = 0,
    HHMPC_DYNMEM_FAIL = 1
} hhmpc_dynmem_error_t;

cJSON *hhmpc_dynmem_get_data(char *fname);


#endif