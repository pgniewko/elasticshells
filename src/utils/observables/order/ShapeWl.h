#ifndef WL_H
#define	WL_H

#include <stdio.h>       /* printf, */
#include <stdlib.h>      /* malloc */
#include <steinhardt.h>  /* */

#include "utils/observables/Observer.h"

class ShapeWl : public Observer
{
    public:
        explicit ShapeWl(const char*, const char*);
        ShapeWl(const ShapeWl& orig);
        virtual ~ShapeWl();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Shell>&, const DomainList&);

    private:
        double calcWl(Shell&);
        static DerivedRegister<ShapeWl> reg;
};

#endif	/* WL_H */

