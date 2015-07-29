#include "include/hhmpctestfunc.h"

/*DEFINITIONEN*/
void change(int a[], int b[]){
    extern int c;
    int i;
    for (i = 0; i < LEN; ++i){
        a[i] = b[i]*b[i];
        printf("%d %d\n", a[i], b[i]);
    }
    c = c*c;
    printf("%d\n", c);
}

int power(int base, int n){
    int i, p;
    p = 1;
    for (i = 1; i <= n; ++i)
        p = p*base;
    return p;
}