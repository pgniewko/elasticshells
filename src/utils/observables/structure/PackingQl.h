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
        double observe(const Box&, const std::vector<Shell>&);

    private:
        double calc_ql(const Box&, const std::vector<Shell>&, unsigned int);
        static DerivedRegister<PackingQl> reg;
};

#endif /* PACKINGQL_H */

