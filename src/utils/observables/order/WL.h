#ifndef WL_H
#define	WL_H

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<steinhardt.h>

#include "Environment.h"
#include "Cell.h"

class WL
{
    public:
        WL();
        WL(const WL& orig);
        virtual ~WL();

        static double calcWl(Cell&, int, double);
        static double calcWl(std::vector<Cell>&, int, double);
    private:

};

#endif	/* WL_H */

