#ifndef HHMPCALG_H
#define HHMPCALG_H

#include "arithmetic.h"

void getBlock(real_t block[], 
              const real_t mtx[], const uint32_t dimension,
              const uint32_t first,const uint32_t size);

void form_Y(void);


#endif