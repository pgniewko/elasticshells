#ifndef TETRAHEDRON_H
#define	TETRAHEDRON_H

#include "Environment.h"
#include "Vector3D.h"

class Tetrahedron
{
    public:
//        Tetrahedron(){}
//        Tetrahedron(Vector3D, Vector3D, Vector3D, Vector3D);
//        Tetrahedron(const Tetrahedron& orig);
//        virtual ~Tetrahedron();
        static double volume(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&);
        static double volume(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&, const double);
        static double volumeSgn(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&);
//        Vector3D a, b, c, d;
//    private:

};

#endif	/* TETRAHEDRON_H */
