#include <stdio.h>
#include <stdlib.h>

#define TEXT    "Zeichen(-kette) eingeben!\n"

int main(/*int argc, char **argv*/) {
    
    int c, d;
    long n[2];
    printf(TEXT);
 
//     c = getchar() == EOF;
    c = getchar();
    printf("%d", c);
    n[0] = 0;
    n[1] = 0;
    while ((d = getchar()) != EOF) {
        if (d == '\n')
            ++n[1];
        putchar(d);
        ++n[0];
//         printf("%d", EOF);
//         printf("%d", c);
//         c = c + 1;
            
    }
    printf("nc%ld nl%ld", n[0], n[1]);
    for (n[0] = 0; getchar() != EOF; ++n[0])
        ;
    printf("%ld", n[0]);
    return 0;
}
