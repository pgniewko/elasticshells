#ifndef TINKER_H
#define	TINKER_H

#include <list>
#include <vector>

#include "Environment.h"
#include "Shell.h"
#include "geometry/Triangle.h"
#include "geometry/Vertex.h"
#include "geometry/VertexTriangle.h"
#include "force/BendingHinge.h"

class Shell;

class Tinker
{
    public:
        Tinker();
        Tinker(const Tinker& orig);
        virtual ~Tinker();
        static void constructVertices(Shell&, std::list<Triangle>&);
        static void constructVTriangles(Shell&, std::list<Triangle>&);
        static void constructTopology(Shell&);
        static void constructBSprings(Shell&);

    private:
        static bool isUnique(std::list<Vector3D>&, Vector3D&, double = constants::epsilon);
        static bool isBSpringUnique(int, int, int, int, Shell&);

        static int getVertex(Shell&, const Vector3D&, double = constants::epsilon);

        static utils::Logger tinker_log;
};

#endif	/* TINKER_H */

