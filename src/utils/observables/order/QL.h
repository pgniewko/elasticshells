#ifndef QL_H
#define	QL_H

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<steinhardt.h>

#include "Environment.h"
#include "Cell.h"
#include "utils/observables/Observer.h"

class QL : public Observer
{
    public:
        QL(const char*, const char*);
        QL(const QL& orig);
        virtual ~QL();

        double observe(Box&, std::vector<Cell>&);
        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);


//        static double calcQl(std::vector<Cell>&, int, double);
    private:
        double calcQl(Cell&);

        static DerivedRegister<QL> reg;
        double l;
        double rc;

};

#endif	/* QL_H */

