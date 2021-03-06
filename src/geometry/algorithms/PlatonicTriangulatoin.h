#ifndef PLATONICTRIANGULATOIN_H
#define	PLATONICTRIANGULATOIN_H

#include <iostream>    /* cout, cin */

#include "Environment.h"
#include "Triangulation.h"

/*
 * This recursive method avoids the common problem of the polar singularity,
 * produced by 2D parameterization methods.
 */
class PlatonicTriangulatoin : public Triangulation
{
    public:
        PlatonicTriangulatoin();
        PlatonicTriangulatoin(int);
        PlatonicTriangulatoin(int, int);
        PlatonicTriangulatoin(const PlatonicTriangulatoin& orig);
        virtual ~PlatonicTriangulatoin();
        std::list<Triangle> triangulate();
        std::list<Triangle> triangulate(double);

    private:
        void create_tetrahedron();
        void create_hexahedron();
        void create_octahedron();
        void create_icosahedron();
        void subdivide();
        int depth;
        int type;

};

#endif	/* PLATONICTRIANGULATOIN_H */