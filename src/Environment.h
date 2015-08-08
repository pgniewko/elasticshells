#ifndef ENVIRONMENT_H
#define	ENVIRONMENT_H

#if defined (_OPENMP)
#include <omp.h>
#endif

#include "random.h"
#include "utils/Logger.h"

#define IMPLIES(x, y) (!(x) || (y))
//#include <assert.h>
//void foo(int array[], int n) {
//  assert(IMPLIES(n > 0, array != NULL));
//  ...

#define STRCMP(a,b) (!strcmp(a,b))
#define SIGN(a) (a >= 0 ? 1 : -1)

#define MAX_CELLS 100
#define MAX_V 400
#define MAX_T 800
#define NEIGH_MAX 20
#define TRIAN_MAX 40
#define NBNEI_MAX 100
#define MAX_IN_DOMAIN 10 // maximum number of particles in a domain
#define MAX_M 50 // maximum number of linked-domains - in every direction
#define MAX_D_NEIGH 27 // MUST BE 27 - ALWAYS !

//TODO: moze warto przerowbic ten plik tak jak tutaj sugeruja:
// http://www.learncpp.com/cpp-tutorial/42-global-variables/

const double E       = 2.71828182845905;
const double PI      = 3.14159265358979;
const double SQRT2   = 1.41421356237310;
const double SQRT3   = 1.73205080756888;
const double SQRT6   = 2.44948974278319;
const double DELTA7  = 0.00000001000000;
const double DELTA14 = 0.00000000000001;
const double EPSILON = 0.01;
const double P3ROOT2 = 1.25992104989487;
const double D4_3    = 1.33333333333333;

const char names[26] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'
                       };

#endif	/* ENVIRONMENT_H */

