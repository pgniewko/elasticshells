#ifndef QL_H
#define	QL_H

#include<stdio.h>
#include<stdlib.h>
#include<steinhardt.h>

#include "utils/observables/Observer.h"

class QL : public Observer
{
    public:
        QL(const char*, const char*);
        QL(const QL& orig);
        virtual ~QL();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);

    private:
        static DerivedRegister<QL> reg;

        double calcQl(Cell&);

};

#endif	/* QL_H */

