#include "Tetrahedron.h"

double Tetrahedron::volume(const Vector3D& a, const Vector3D& b, const Vector3D& c, const Vector3D& cm)
{
    return fabs((a - cm) * cross(b - cm, c - cm)) / 6.0;
}

double Tetrahedron::volume(const Vector3D& a, const Vector3D& b, const Vector3D& c, const Vector3D& cm, const double eps)
{
    if (eps == 0)
    {
        return volume(a, b, c, cm);
    }

    Vector3D AD = a - cm;
    Vector3D BD = b - cm;
    Vector3D CD = c - cm;
    double nAD = AD.length() + eps;
    double nBD = BD.length() + eps;
    double nCD = CD.length() + eps;

    AD.set_length(nAD);
    BD.set_length(nBD);
    CD.set_length(nCD);

    return ( fabs(AD * cross(BD, CD)) / 6.0 );
}

double Tetrahedron::volume_sgn(const Vector3D& a, const Vector3D& b, const Vector3D& c, const Vector3D& cm)
{
    double volume = (a - cm) * cross((b - cm), (c - cm));
    return  SIGN( volume );
}