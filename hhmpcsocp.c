#include "hhmpcsocp.h"


static void hhmpc_copy_data(struct hhmpc_term *dest, struct hhmpc_term *src);

void hhmpc_socp_form_problem(struct hhmpc_socp *socp)
{
    /* Needed if parametric Matrices are used*/
}

void hhmpc_copy_data(struct hhmpc_term *dest, struct hhmpc_term *src)  {
    int j;
    for (j=0; j<(dest->cols*dest->rows); j++) {
        dest->data[j] = src->data[j];
    }
    return;
}

