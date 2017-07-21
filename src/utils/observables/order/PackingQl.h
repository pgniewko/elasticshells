#ifndef PACKINGQL_H
#define PACKINGQL_H

#include <stdio.h>       /* printf, */
#include <stdlib.h>      /* malloc */
#include <steinhardt.h>  /* */

#include "utils/observables/Observer.h"

class PackingQl : public Observer
{
    public:
        explicit PackingQl(const char*, const char*);
        PackingQl(const PackingQl& orig);
        virtual ~PackingQl();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&, const DomainList&);

    private:
        double calcQl(const Box&, std::vector<Cell>&, unsigned int, const DomainList&);
        static DerivedRegister<PackingQl> reg;
};

#endif /* PACKINGQL_H */

