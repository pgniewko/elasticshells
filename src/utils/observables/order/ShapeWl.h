#ifndef WL_H
#define	WL_H

#include<stdio.h>       /* printf, */
#include<stdlib.h>      /* malloc */
#include<steinhardt.h>  /* */

#include "utils/observables/Observer.h"

class ShapeWl : public Observer
{
    public:
        explicit ShapeWl(const char*, const char*);
        ShapeWl(const ShapeWl& orig);
        virtual ~ShapeWl();

        double observe(const Box&, std::vector<Cell>&);
        void set_params(const int, std::vector<std::string>);

    private:
        double calcWl(Cell&);
        static DerivedRegister<ShapeWl> reg;
};

#endif	/* WL_H */

