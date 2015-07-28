#ifndef WL_H
#define	WL_H

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<steinhardt.h>

#include "utils/observables/Observer.h"

class WL : public Observer
{
    public:
        WL(const char*, const char*);
        WL(const WL& orig);
        virtual ~WL();

        double observe(Box&, std::vector<Cell>&);
        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);
        
    private:
        double calcWl(Cell&);
        static DerivedRegister<WL> reg;
};

#endif	/* WL_H */

