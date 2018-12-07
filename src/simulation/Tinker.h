#ifndef TINKER_H
#define	TINKER_H

#include <list>
#include <vector>

#include "Environment.h"
#include "Shell.h"
#include "geometry/Triangle.h"
#include "geometry/Vertex.h"
#include "geometry/Element.h"
#include "geometry/Hinge.h"

class Shell;

class Tinker
{
    public:
        Tinker();
        Tinker(const Tinker& orig);
        virtual ~Tinker();
        static void construct_vertices(Shell&, std::list<Triangle>&);
        static void construct_elements(Shell&, std::list<Triangle>&);
        static void construct_topology(Shell&);
        static void construct_hinges(Shell&);

    private:
        static bool is_unique(std::list<Vector3D>&, Vector3D&, double = constants::epsilon);
        static bool is_hinge_unique(int, int, int, int, Shell&);

        static int get_vertex(Shell&, const Vector3D&, double = constants::epsilon);

        static utils::Logger tinker_log;
};

#endif	/* TINKER_H */

