#ifndef HHMPCUSEFULL_H
#define HHMPCUSEFULL_H

#include <stdio.h>
#include "arithmetic.h"

void eye(real_t mtx[], const uint32_t dimension);

void zeroes(real_t mtx[], const uint32_t length);

void print_mtx(const real_t mtx[], const uint32_t rows, const uint32_t colums);


#endif