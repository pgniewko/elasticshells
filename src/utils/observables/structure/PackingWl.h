#ifndef PACKINGWL_H
#define PACKINGWL_H

#include <stdio.h>       /* printf, */
#include <stdlib.h>      /* malloc */
#include <steinhardt.h>  /* */

#include "utils/observables/Observer.h"

class PackingWl : public Observer
{
    public:
        explicit PackingWl(const char*, const char*);
        PackingWl(const PackingWl& orig);
        virtual ~PackingWl();

        double observe(const Box&, const std::vector<Shell>&);

    private:
        double calc_wl(const Box&, const std::vector<Shell>&, unsigned int);
        static DerivedRegister<PackingWl> reg;
};

#endif /* PACKINGWL_H */

