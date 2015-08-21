#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */
#include "include/cjson.h"

#include "include/hhmpcsocpdynmem.h"

/* Extern function definitions */

struct hhmpc_socp *hhmpc_socp_allocate_former(void)
{
    
    struct hhmpc_socp *socp = (struct hhmpc_socp*)malloc(sizeof(struct hhmpc_socp));
    
    return socp;
}