#ifndef PACKINGWL_H
#define PACKINGWL_H

#include<stdio.h>       /* printf, */
#include<stdlib.h>      /* malloc */
#include<steinhardt.h>  /* */

#include "utils/observables/Observer.h"

class PackingWl  : public Observer
{
    public:
        explicit PackingWl(const char*, const char*);
        PackingWl(const PackingWl& orig);
        virtual ~PackingWl();
    
        double observe(const Box&, std::vector<Cell>&);
        void set_params(const int, std::vector<std::string>);
        
    private:
        static DerivedRegister<PackingWl> reg;

};

#endif /* PACKINGWL_H */

