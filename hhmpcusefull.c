#include "include/hhmpcusefull.h"


void eye(real_t mtx[], const uint32_t dim)
{
    uint32_t i, j;
    for (i = 0; i < dim; i++){
        for (j = 0; j < dim; j++)
            mtx[i*dim+j] = (i == j) ? 1. : 0.;
    }
}

void zeroes(real_t mtx[], const uint32_t l)
{
    uint32_t i;
    for (i = 0; i < l; i++)
        mtx[i] = 0.;
}

void print_mtx(const real_t mtx[], const uint32_t rows, const uint32_t cols)
{
//     uint32_t i, j;
// //     printf("[\n");
//     for (i = 0; i < rows; i++){
//         for (j = 0; j < cols; j++)
//             printf("\t%.3f,", mtx[i*cols + j]);
// //         printf("\n");
//     }
//     printf("]\n");
}
