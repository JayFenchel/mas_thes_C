#include "include/hhmpctestfunc.h"

#define ALLOCSIZE 10

static float allocbuf[ALLOCSIZE];
static float *allocp = allocbuf; /*nächste freie Position, zeigt zum Beginn auf den Anfang von allocbuf*/

/*DEFINITIONEN*/

float *alloc(int n){ /*liefert Zeiger auf Platz für n Zeichen*/
    if (allocbuf + ALLOCSIZE - allocp >= n) {
        allocp += n;
        return allocp - n;
    } else        
    return 0;
}

void afree(float *p){
    if (p >= allocbuf && p < allocbuf + ALLOCSIZE)
        allocp = p;
}

void aprint(void){
    int i;
    for (i=0; i<ALLOCSIZE; ++i){
        printf("%f\n", allocbuf[i]);
    }
}


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

int sum(int a, int b) {

    return a + b;
}