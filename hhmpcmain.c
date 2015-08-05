#include <stdio.h>
#include <stdlib.h>
#include "hhmpctestfunc.h"

#define PPP     printf("help\n")


int main(void){
    extern int c;
    int i, *p1, *p2, v[LEN];
    float *pc;
    
    pc = alloc(6);
    
    pc[0] = 3.4;
    *(pc+5) = 2.5;
    pc[6] = 7.;
    pc = alloc(2);
    pc[0] = 4.4;
    afree(pc-4);
    
    pc = alloc(2);
    pc[0] = 4.3;
    
    aprint();
    
    i = 5;
    p1 = &i;
    printf("%d \n", *p1);
    
    for (i = 0; i < LEN; ++i){
        printf("%d %d %d\n", i, power(2, i), power(-3, i));
        v[i] = power(2, i);
    }
    
    p2 = v;
    for (i=0; i < LEN; ++i){
        printf("%d ", v[i]);
        printf("%d ", p2[i]);
    }
    printf("\n");
    for (i=0; i < LEN; ++i){
        printf("%d ", *v);
        printf("%d ", *(v+i));
    }
    printf("\n");
    for (i=0; i < LEN; ++i){
        printf("%d ", *p2);
        printf("%d ", *(p2+i));
    }
    printf("\n");
    for (i=0; i < LEN; ++i){
        printf("%d ", *(p2++));
    }
    printf("\n");
    
    c = SQUARE(1+1+1);
    change(v, v);
    for (i = 0; i < LEN; ++i)
        printf("%d\n", c);
    PPP;
    return 0;
}


// int main(/*int argc, char **argv*/) {
//     
//     int c, d;
//     long n[2];
//     printf(TEXT\n);
//  
// //     c = getchar() == EOF;
//     c = getchar();
//     printf("%d", c);
//     n[0] = 0;
//     n[1] = 0;
//     while ((d = getchar()) != EOF) {
//         if (d == '\n')
//             ++n[1];
//         putchar(d);
//         ++n[0];
// //         printf("%d", EOF);
// //         printf("%d", c);
// //         c = c + 1;
//             
//     }
//     printf("nc%ld nl%ld", n[0], n[1]);
//     for (n[0] = 0; getchar() != EOF; ++n[0])
//         ;
//     printf("%ld", n[0]);
//     return 0;
// }

