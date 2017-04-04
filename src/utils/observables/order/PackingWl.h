#ifndef PACKINGWL_H
#define PACKINGWL_H

#include <stdio.h>       /* printf, */
#include <stdlib.h>      /* malloc */
#include <steinhardt.h>  /* */

#include "utils/observables/Observer.h"

class PackingWl  : public Observer
{
    public:
        explicit PackingWl(const char*, const char*);
        PackingWl(const PackingWl& orig);
        virtual ~PackingWl();
    
        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&, const DomainList&);
        
    private:
        double calcWl(const Box&, std::vector<Cell>&, unsigned int, const DomainList&);
        static DerivedRegister<PackingWl> reg;
};

#endif /* PACKINGWL_H */

