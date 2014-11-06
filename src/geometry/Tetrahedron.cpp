#include "Tetrahedron.h"

Tetrahedron::Tetrahedron(Vector3D m, Vector3D n, Vector3D o, Vector3D p) : a(m), b(n), c(o), d(p) {}

Tetrahedron::Tetrahedron(const Tetrahedron& orig) : a(orig.a), b(orig.b), c(orig.c), d(orig.d) {}

Tetrahedron::~Tetrahedron()
{
}

double Tetrahedron::volume() const
{
    Vector3D AD = a - d;
    Vector3D BD = b - d;
    Vector3D CD = c - d;
    Vector3D CcrossD = cross(BD, CD);
    double volume = fabs(AD * CcrossD);
    volume /= 6.0;
    return volume;
}

double Tetrahedron::volumeSgn() const
{
    Vector3D AD = a - d;
    Vector3D BD = b - d;
    Vector3D CD = c - d;
    Vector3D CcrossD = cross(BD, CD);
    double volume = AD * CcrossD;
    //volume /= 6.0;
    return  SIGN( volume );
}