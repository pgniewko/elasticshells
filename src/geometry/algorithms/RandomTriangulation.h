#ifndef RANDOMTRIANGULATION_H
#define RANDOMTRIANGULATION_H

#include <rndmesh.h>
#include <iostream>    /* cout, cin */
#include <math.h>

#include "Environment.h"
#include "Triangulation.h"


class RandomTriangulation : public Triangulation
{
    public:
        RandomTriangulation(int, int, double, double, double);
        RandomTriangulation(const RandomTriangulation& orig);
        virtual ~RandomTriangulation();

        std::list<Triangle> triangulate();
        std::list<Triangle> triangulate(double);
        std::list<Triangle> triangulate(double, double, double, int);

    private:

        int n_steps;
        int n_anneals;
        double T_min;
        double T_max;
        double r_vertex;

};

#endif /* RANDOMTRIANGULATION_H */

