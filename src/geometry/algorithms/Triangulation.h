#ifndef TRIANGULATION_H
#define	TRIANGULATION_H

#include <list>         /* list*/
#include <iostream>

#include "utils/utils.h"
#include "Environment.h"
#include "geometry/Triangle.h"

class Triangulation
{
    public:
        Triangulation();
        Triangulation(const Triangulation& orig);
        virtual ~Triangulation();
        virtual std::list<Triangle> triangulate() = 0;
    protected:
        std::list<Triangle> tris;

};

#endif	/* TRIANGULATION_H */
