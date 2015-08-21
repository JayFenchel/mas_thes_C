#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */
#include "include/cjson.h"

#include "include/hhmpcsocpdynmem.h"

/* Extern function definitions */

struct hhmpc_socp *hhmpc_socp_allocate_former(void)
{
    int i;
    struct hhmpc_term *t;
    struct hhmpc_socp *socp = (struct hhmpc_socp*)malloc(sizeof(struct hhmpc_socp));
    if (NULL == socp) {return NULL;}
    
    /* constants */
    t = (struct hhmpc_term*)calloc(HHMPC_PAR_NUM, sizeof(struct hhmpc_term));
    if (NULL == t) {return NULL;}
    for (i=0; i<HHMPC_PAR_NUM; i++) {
        socp->constant[i] = &t[i];
    }
    
    /* the evaluated problem itself */
    socp->prb = (struct hhmpc_socp_prb*)malloc(sizeof(struct hhmpc_socp_prb));
    if (NULL == socp->prb) {return NULL;}
        socp->prb->q = socp->constant[HHMPC_Q_KL];
        socp->prb->r = socp->constant[HHMPC_R_KL];
    
    return socp;
}