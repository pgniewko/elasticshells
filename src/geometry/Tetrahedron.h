#ifndef TETRAHEDRON_H
#define	TETRAHEDRON_H

#include "Environment.h"
#include "Vector3D.h"

class Tetrahedron
{
    public:
        static double volume(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&);
        static double volume(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&, const double);
        static double volumeSgn(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&);
};

#endif	/* TETRAHEDRON_H */
