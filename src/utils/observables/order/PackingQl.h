#ifndef PACKINGQL_H
#define PACKINGQL_H

#include<stdio.h>       /* printf, */
#include<stdlib.h>      /* malloc */
#include<steinhardt.h>  /* */

#include "utils/observables/Observer.h"

class PackingQl : public Observer
{
    public:
        explicit PackingQl(const char*, const char*);
        PackingQl(const PackingQl& orig);
        virtual ~PackingQl();
        
        double observe(const Box&, std::vector<Cell>&);
        void set_params(const int, std::vector<std::string>);
        
private:
    
    static DerivedRegister<PackingQl> reg;

};

#endif /* PACKINGQL_H */

