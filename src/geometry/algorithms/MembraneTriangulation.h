#ifndef MEMBRANETRIANGULATION_H
#define	MEMBRANETRIANGULATION_H

#include <iostream>    /* cout, cin */

#include "Environment.h"
#include "Triangulation.h"

class MembraneTriangulation : public Triangulation
{
    public:
        MembraneTriangulation();
        MembraneTriangulation(const MembraneTriangulation& orig);
        virtual ~MembraneTriangulation();
        std::list<Triangle> triangulate();
        std::list<Triangle> triangulate(double, double, int);
    private:
        std::list<Triangle> hexagon;
        void createHex(double, double, double);
        void subdivide();


};

#endif	/* MEMBRANETRIANGULATION_H */

