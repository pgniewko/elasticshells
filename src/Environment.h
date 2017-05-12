#ifndef ENVIRONMENT_H
#define	ENVIRONMENT_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <string.h>
#include <limits>
#include <float.h>
#include "random.h"
#include "fastmath.h"
#include "utils/Logger.h"

#define IMPLIES(x, y) (!(x) || (y))

#define STRCMP(a,b) (!strcmp(a,b))
#define SIGN(a) (a >= 0 ? 1 : -1)

#define MAX_CELLS 100
#define MAX_V 2501
#define MAX_T 5000
#define NEIGH_MAX 12
#define TRIAN_MAX 24
#define MAX_M 50 // maximum number of linked-domains - in each direction
#define MAX_D_NEIGH 26 // MUST BE 26 - ALWAYS !

namespace constants
{
    const double e       = 2.71828182845905;
    const double pi      = 3.14159265358979;
    const double sqrt2   = 1.41421356237310;
    const double sqrt3   = 1.73205080756888;
    const double sqrt6   = 2.44948974278319;
    const double delta7  = 0.00000001;
    const double delta14 = 0.00000000000001;
    const double epsilon = 0.001;
    const double p3root2 = 1.25992104989487;
    const double d4_3    = 1.33333333333333;
    const double d2_5    = 0.40000000000000;
    const char names[26] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                            'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'
                           };
}

namespace analyze
{
        double MIN_FORCE_SQ = 0.0;
        double FORCE_FRAC = 0.0;
}

#endif	/* ENVIRONMENT_H */

