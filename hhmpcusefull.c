#include "include/hhmpcusefull.h"


void zeroes(real_t mtx[], const uint32_t l)
{
    uint32_t i;
    for (i = 0; i < l; i++)
        mtx[i] = 0;
}

void print_mtx(const real_t mtx[], const uint32_t rows, const uint32_t cols)
{
    uint32_t i, j;
    printf("[\n");
    for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++)
            printf("\t%.4f", mtx[i*cols + j]);
        printf("\n");
    }
    printf("]\n");
}
