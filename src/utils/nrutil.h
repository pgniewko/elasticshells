#ifndef NRUTIL_H
#define	NRUTIL_H

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#define SIGN2(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

//static double minarg1, minarg2;
//#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
//        (minarg1) : (minarg2))

void nrerror(char error_text[]);
double* darray(long n);
void free_darray(double* da);

#endif	/* NRUTIL_H */
