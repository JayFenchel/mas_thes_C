#ifndef MPCINC_FGMDYNMEM_H
#define MPCINC_FGMDYNMEM_H

#include "mpcincdynmem.h"
#include "mpcinccvp.h"
#include "mpcincfgm.h"

extern struct mpcinc_fgm *mpcinc_fgm_allocate_solver(void);

extern mpcinc_dynmem_error_t mpcinc_fgm_setup_solver(
                struct mpcinc_fgm *fgm,
                struct mpcinc_cvp_prb *prb, char *fname);

#endif /* MPCINC_FGMDYNMEM_H */
