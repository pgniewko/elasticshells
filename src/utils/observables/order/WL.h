#ifndef WL_H
#define	WL_H

#include<stdio.h>       /* printf, */
#include<stdlib.h>      /* malloc */
#include<steinhardt.h>  /* */

#include "utils/observables/Observer.h"

class WL : public Observer
{
    public:
        explicit WL(const char*, const char*);
        WL(const WL& orig);
        virtual ~WL();

        double observe(const Box&, std::vector<Cell>&);
        void set_params(const int, std::vector<std::string>);

    private:
        double calcWl(Cell&);
        static DerivedRegister<WL> reg;
};

#endif	/* WL_H */

