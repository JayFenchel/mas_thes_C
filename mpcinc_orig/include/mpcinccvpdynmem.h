#ifndef MPCINC_CVPDYNMEM_H
#define MPCINC_CVPDYNMEM_H

#include "mpcincdynmem.h"
#include "mpcinccvp.h"

extern mpcinc_dynmem_error_t mpcinc_cvp_setup_former(struct mpcinc_cvp *cvp, char *fname);
extern struct mpcinc_cvp *mpcinc_cvp_allocate_former(void);

#endif /* MPCINC_CVPDYNMEM_H */
