#include <stdlib.h>  /* malloc */
#include "include/hhmpcipmdynmem.h"

/* Extern function definitions */

struct hhmpc_ipm *hhmpc_ipm_allocate_solver(void)
{
    struct hhmpc_ipm *ipm = (struct hhmpc_ipm*)malloc(sizeof(struct hhmpc_ipm));
    if (NULL == ipm) {return NULL;}
    struct hhmpc_ipm_conf *conf = (struct hhmpc_ipm_conf*)malloc(sizeof(struct hhmpc_ipm_conf));
    if (NULL == conf) {return NULL;}
    ipm->conf = conf;
    
    return ipm;
}
