#ifndef QL_H
#define	QL_H

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<steinhardt.h>

#include "Environment.h"
#include "Cell.h"

class QL {
public:
    QL();
    QL(const QL& orig);
    virtual ~QL();
    
    static double calcQl(Cell, int, double);
private:

};

#endif	/* QL_H */

