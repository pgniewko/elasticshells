#ifndef QL_H
#define	QL_H

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<steinhardt.h>

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
        
    private:
        static DerivedRegister<QL> reg;
        
        double calcQl(Cell&);
        double l;
        double rc;

};

#endif	/* QL_H */

