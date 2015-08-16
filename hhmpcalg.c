#include "include/hhmpcalg.h"


void form_Y(void)
{
    
}

void getBlock(real_t blk[],
              const real_t mtx[], const uint32_t dim,
              const uint32_t fst, const uint32_t s)
{
    uint32_t i, j; /* loop counters */
    
    for (i = 0; i < s; i++){
        for (j = 0; j < s; j++){
            blk[i*s+j] = mtx[(i+fst)*dim+j+fst];
        }
    }
}