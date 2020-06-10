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
        std::list<Triangle> diamond1;
        std::list<Triangle> diamond2;
        void create_hex(double, double);
        void put_two_triangles(double);
        void subdivide();
};

#endif	/* MEMBRANETRIANGULATION_H */

