#ifndef LU_H_
#define LU_H_

#include <math.h>
#include <iostream>
using namespace std;

// prototypes
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);

#endif /*LU_H_*/
