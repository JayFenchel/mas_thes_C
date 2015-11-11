#ifndef STATIC_DATA_H
#define STATIC_DATA_H

/* Define desired example here. */
// #define HHMPC_SOCPCONDTEST  /* enable for HHMPC_QPCONDTESTs, too */
// #define HHMPC_SOCPCONDTEST5
// #define HHMPC_QPCONDTEST5
// #define HHMPC_QPCONDTEST10
// #define HHMPC_QPCONDTEST20
// #define HHMPC_QPCONDTEST30
// #define HHMPC_QPTEST
#define HHMPC_QPSMALLTEST
// #define HHMPC_QPSMALLTEST5
// #define HHMPC_QPSMALLTEST10
// #define HHMPC_QPSMALLTEST20
#define HHMPC_QPSMALLTEST30
// #define HHMPC_SOCPTEST
// #define HHMPC_SOCPHARDTEST
// #define HHMPC_SOCPSOFTTEST

#include "hhmpcsocp.h"
#include "mpccvp.h"
#include "hhmpcipm.h"
#include "mpcincmtxops.h"
#include "hhmpcusefull.h"

#define HHMPC_SIMPOINTS 60

#ifdef HHMPC_SOCPCONDTEST
#define HHMPC_NB_ROWSFU 4
#define HHMPC_NB_ROWSFF 2

#ifdef HHMPC_SOCPCONDTEST5

#define HHMPC_HORIZON 5
#define HHMPC_SV 0
#define HHMPC_CV 1
#define HHMPC_OV 1
#define HHMPC_OS 5
#define HHMPC_DS 0
#define HHMPC_NB_LCONSTR 10
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 5
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 0

#endif  /* HHMPC_SOCPCONDTEST5 */

#ifdef HHMPC_QPCONDTEST5

#define HHMPC_HORIZON 5
#define HHMPC_SV 0
#define HHMPC_CV 1
#define HHMPC_OV 1
#define HHMPC_OS 5
#define HHMPC_DS 0
#define HHMPC_NB_LCONSTR 20
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 0
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 0

#endif  /* HHMPC_QPCONDTEST5 */

#ifdef HHMPC_QPCONDTEST10

#define HHMPC_HORIZON 10
#define HHMPC_SV 0
#define HHMPC_CV 1
#define HHMPC_OV 1
#define HHMPC_OS 10
#define HHMPC_DS 0
#define HHMPC_NB_LCONSTR 40
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 0
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 0

#endif  /* HHMPC_QPCONDTEST10 */

#ifdef HHMPC_QPCONDTEST20

#define HHMPC_HORIZON 20
#define HHMPC_SV 0
#define HHMPC_CV 1
#define HHMPC_OV 1
#define HHMPC_OS 20
#define HHMPC_DS 0
#define HHMPC_NB_LCONSTR 80
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 0
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 0

#endif  /* HHMPC_QPCONDTEST20 */

#ifdef HHMPC_QPCONDTEST30

#define HHMPC_HORIZON 30
#define HHMPC_SV 0
#define HHMPC_CV 1
#define HHMPC_OV 1
#define HHMPC_OS 30
#define HHMPC_DS 0
#define HHMPC_NB_LCONSTR 120
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 0
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 0

#endif  /* HHMPC_QPCONDTEST30 */

#endif

#ifdef HHMPC_QPTEST

#define HHMPC_HORIZON 5
#define HHMPC_SV 30
#define HHMPC_CV 1
#define HHMPC_OV 31
#define HHMPC_OS 155
#define HHMPC_DS 150
#define HHMPC_NB_LCONSTR 22
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 0
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 0

#endif

#ifdef HHMPC_QPSMALLTEST5

#define HHMPC_HORIZON 5
#define HHMPC_SV 5
#define HHMPC_CV 1
#define HHMPC_OV 6
#define HHMPC_OS 30
#define HHMPC_DS 25
#define HHMPC_NB_LCONSTR 22
#define HHMPC_NB_ROWSFU 4
#define HHMPC_NB_ROWSFF 2
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 0
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 0

#endif

#ifdef HHMPC_QPSMALLTEST10

#define HHMPC_HORIZON 10
#define HHMPC_SV 5
#define HHMPC_CV 1
#define HHMPC_OV 6
#define HHMPC_OS 60
#define HHMPC_DS 50
#define HHMPC_NB_LCONSTR 42
#define HHMPC_NB_ROWSFU 4
#define HHMPC_NB_ROWSFF 2
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 0
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 0

#endif

#ifdef HHMPC_QPSMALLTEST20

#define HHMPC_HORIZON 20
#define HHMPC_SV 5
#define HHMPC_CV 1
#define HHMPC_OV 6
#define HHMPC_OS 120
#define HHMPC_DS 100
#define HHMPC_NB_LCONSTR 82
#define HHMPC_NB_ROWSFU 4
#define HHMPC_NB_ROWSFF 2
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 0
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 0

#endif

#ifdef HHMPC_QPSMALLTEST30

#define HHMPC_HORIZON 30
#define HHMPC_SV 5
#define HHMPC_CV 1
#define HHMPC_OV 6
#define HHMPC_OS 180
#define HHMPC_DS 150
#define HHMPC_NB_LCONSTR 122
#define HHMPC_NB_ROWSFU 4
#define HHMPC_NB_ROWSFF 2
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 0
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 0

#endif

#ifdef HHMPC_SOCPTEST

#define HHMPC_HORIZON 5
#define HHMPC_SV 30
#define HHMPC_CV 1
#define HHMPC_OV 31
#define HHMPC_OS 155
#define HHMPC_DS 150
#define HHMPC_NB_LCONSTR 10  /* Only for input constraints */
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 5  /* Here SOCCs are used */
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 0

#endif

#ifdef HHMPC_SOCPHARDTEST

#define HHMPC_HORIZON 5
#define HHMPC_SV 30
#define HHMPC_CV 1
#define HHMPC_OV 31
#define HHMPC_OS 155
#define HHMPC_DS 150
#define HHMPC_NB_LCONSTR 22  /* Only for input constraints */
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 5  /* Here SOCCs are used */
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 0

#endif

#ifdef HHMPC_SOCPSOFTTEST

#define HHMPC_HORIZON 5
#define HHMPC_SV 30
#define HHMPC_CV 1
#define HHMPC_OV 31
#define HHMPC_OS 155
#define HHMPC_DS 150
#define HHMPC_NB_LCONSTR 10  /* Only for input constraints */
#define HHMPC_NB_QC 0
#define HHMPC_NB_SOCC 0  /* Here SOCCs are used */
#define HHMPC_NB_IEQ (HHMPC_NB_LCONSTR+HHMPC_NB_QC+HHMPC_NB_SOCC)
#define HHMPC_NB_SOFT 22

#endif

extern void form_ipm(struct hhmpc_ipm *ipm, struct hhmpc_socp_prb *prb);

extern void form_socp_H(struct hhmpc_socp *socp);

extern void form_socp_with_cvp(struct hhmpc_socp *socp, struct mpc_cvp_prb *cvp_prb);

extern void form_socp(struct hhmpc_socp *socp);


#endif /* STATIC_DATA_H */

