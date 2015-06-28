#include "OsmoticForce.h"

OsmoticForce::OsmoticForce() {}

OsmoticForce::OsmoticForce(const OsmoticForce& orig) {}

OsmoticForce::~OsmoticForce() {}

Vector3D OsmoticForce::calcForce(const Vector3D& va, const Vector3D& vb, const Vector3D& vc, const Vector3D& vd, double nRT, double vol, const double dp)
{
    Vector3D BD = vb - vd;
    Vector3D CD = vc - vd;
    Vector3D f = cross(BD, CD) / 6;

    if (volumeFlag)
    {
        f *= ( nRT / vol );
    }
    else
    {
        f *= dp;
    }

    Tetrahedron tetra(va, vb, vc, vd);
    return tetra.volumeSgn() * f;
}
void OsmoticForce::setVolumeFlag(bool flag)
{
    volumeFlag = flag;
}

void OsmoticForce::setEpsilon(double eps)
{
    epsilon = eps;
}

double OsmoticForce::getEpsilon()
{
    return epsilon;
}

const bool OsmoticForce::getFlag()
{
    return volumeFlag;
}

//double OsmoticForce::epsilon = 0.0;
double OsmoticForce::epsilon = 0.4;
bool OsmoticForce::volumeFlag = false;