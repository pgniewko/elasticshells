#ifndef TETRAHEDRON_H
#define	TETRAHEDRON_H

#include "Environment.h"
#include "Vector3D.h"

class Tetrahedron
{
    public:
//    Tetrahedron();
        Tetrahedron(Vector3D, Vector3D, Vector3D, Vector3D);
        Tetrahedron(const Tetrahedron& orig);
        virtual ~Tetrahedron();
        double volume() const;
        double volume(double) const;
        double volumeSgn() const;
        Vector3D a, b, c, d;
    private:

};

#endif	/* TETRAHEDRON_H */
