#ifndef SHAPEQL_H
#define	SHAPEQL_H

#include<stdio.h>
#include<stdlib.h>
#include<steinhardt.h>

#include "utils/observables/Observer.h"

class ShapeQl : public Observer
{
    public:
        explicit ShapeQl(const char*, const char*);
        ShapeQl(const ShapeQl& orig);
        virtual ~ShapeQl();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);

    private:
        double calcQl(Cell&);
        static DerivedRegister<ShapeQl> reg;
};

#endif	/* SHAPEQL_H */

