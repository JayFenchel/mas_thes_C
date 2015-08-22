#ifndef HHMPCIPMDYNMEM_H
#define HHMPCIPMDYNMEM_H


#include "hhmpcdynmem.h"
#include "hhmpcsocp.h"
#include "hhmpcipm.h"

extern struct hhmpc_ipm *hhmpc_ipm_allocate_solver(void);

extern hhmpc_dynmem_error_t hhmpc_ipm_setup_solver(struct hhmpc_ipm *ipm,
                                                   struct hhmpc_socp_prb *prb,
                                                   char *fname);


#endif /* HHMPCIPMDYNMEM_H */