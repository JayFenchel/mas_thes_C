#ifndef STATIC_DATA_H
#define STATIC_DATA_H


#include "hhmpcsocp.h"
#include "hhmpcipm.h"
#include "mpcincmtxops.h"
#include "hhmpcusefull.h"


#define HHMPC_HORIZON 5
#define HHMPC_SV 30
#define HHMPC_CV 1
#define HHMPC_OV 31
#define HHMPC_OS 155
#define HHMPC_DS 150
#define HHMPC_NB_LCONSTR 34
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 5
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 34

extern void form_ipm(struct hhmpc_ipm *ipm, struct hhmpc_socp_prb *prb);

extern void form_socp(struct hhmpc_socp *socp);


#endif /* STATIC_DATA_H */

