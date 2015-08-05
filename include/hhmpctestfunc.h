#ifndef HHMPCTESTFUNC_H
#define HHMPCTESTFUNC_H

#include <stdio.h>

#define LEN     10
#define TEXT    "Zeichen(-kette) eingeben!"
#define SQUARE(x)       (x)*(x)

int c;

/*DEKLARATIONEN*/
float *alloc(int n_speicher_zellen);

void afree(float *buff_auf_p_setzten);

void aprint(void);

void change(int wrong[], int right[]);

int power(int m, int n);

int sum(int a, int b);


#endif