#ifndef HHMPCSOCPDYNMEM_H
#define HHMPCSOCPDYNMEM_H


#include "hhmpcdynmem.h"
#include "hhmpcsocp.h"

extern hhmpc_dynmem_error_t hhmpc_socp_setup_former(struct hhmpc_socp *socp,
                                                    char *fname);
extern struct hhmpc_socp *hhmpc_socp_allocate_former(void);


#endif /* HHMPCSOCPDYNMEM_H */