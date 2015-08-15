#include "include/hhmpcmath.h"

real_t mtx_out(real_t mtx[], uint32_t rows, uint32_t cols, real_t vec[]){
    real_t out[2];
    mpcinc_mtx_multiply_mtx_vec(out, mtx, vec, rows, cols);
    return out[0];
}

int simple_sum(int a, int b) {

    return a + b;
}