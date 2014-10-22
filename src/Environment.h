#ifndef ENVIRONMENT_H
#define	ENVIRONMENT_H

#include "random.h"
#include "utils/Logger.h"

#define STRCMP(a,b) (!strcmp(a,b))

#define MAX_CELLS 100
#define MAX_V 200
#define MAX_T 500
#define NEIGH_MAX 10
#define TRIAN_MAX 10
#define NBNEI_MAX 100
#define MAX_IN_DOMAIN 20 // maximum number of particles in a domain
#define MAX_M 20 // maximum number of linked-domains

const double E       = 2.71828182845905;
const double PI      = 3.14159265358979;
const double SQRT2   = 1.41421356237310;
const double SQRT3   = 1.73205080756888;
const double SQRT6   = 2.44948974278319;
const double DELTA7  = 0.00000001000000;
const double DELTA14 = 0.00000000000001;
const double EPSILON = 0.25;
const double P3ROOT2 = 1.25992104989487;

const char names[26] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};

#endif	/* ENVIRONMENT_H */

