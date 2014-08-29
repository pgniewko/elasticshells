#include "OsmoticForce.h"

OsmoticForce::OsmoticForce() {}

OsmoticForce::OsmoticForce(const OsmoticForce& orig) {}

OsmoticForce::~OsmoticForce() {}

static Vector3D OsmoticForce::calcForce(const Vector3D& va, const Vector3D& vb, const Vector3D& vc, const Vector3D& vd, const double dp)
{
        Vector3D BD = vb - vd;
        Vector3D CD = vc - vd;
        Vector3D f = -dp * cross(BD, CD) / 6;
        return f;
}