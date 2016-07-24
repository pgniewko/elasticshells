#ifndef NRUTIL_H
#define	NRUTIL_H

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#define SIGN2(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define ARRAY_SIZE(array) ( sizeof(array) / sizeof(array[0]) )
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?(maxarg1) : (maxarg2))

void nrerror(char error_text[]);
double* darray(long n);
void free_darray(double* da);

#endif	/* NRUTIL_H */
