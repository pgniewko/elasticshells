#include "Tetrahedron.h"

double Tetrahedron::volume(const Vector3D& a, const Vector3D& b, const Vector3D& c, const Vector3D& d)
{
    return fabs((a-d) * cross(b-d, c-d)) / 6.0;
}

double Tetrahedron::volume(const Vector3D& a, const Vector3D& b, const Vector3D& c, const Vector3D& d, const double eps)
{
    Vector3D AD = a - d;
    Vector3D BD = b - d;
    Vector3D CD = c - d;
    double nAD = AD.length() + eps;
    double nBD = BD.length() + eps;
    double nCD = CD.length() + eps;
    AD.set_length(nAD);
    BD.set_length(nBD);
    CD.set_length(nCD);
    Vector3D CcrossD = cross(BD, CD);
    double volume = fabs(AD * CcrossD);
    volume /= 6.0;
    return volume;
}

double Tetrahedron::volumeSgn(const Vector3D& a, const Vector3D& b, const Vector3D& c, const Vector3D& d)
{
    double volume = (a - d) * cross((b-d), (c-d));
    return  SIGN( volume );   
}