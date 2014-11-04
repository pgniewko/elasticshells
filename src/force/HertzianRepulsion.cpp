#include "HertzianRepulsion.h"

HertzianRepulsion::HertzianRepulsion() {}

HertzianRepulsion::HertzianRepulsion(const HertzianRepulsion& orig) {}

HertzianRepulsion::~HertzianRepulsion() {}

Vector3D HertzianRepulsion::calcForce(Vector3D& dij, double R, double E)
{
    double eps = R - dij.length();
    if (eps > 0)
    {
        double fmagn = D4_3 * E * pow(R, 0.5) * pow(eps, 1.5);
        Vector3D f = fmagn * (dij / dij.length()) ;
        return f;
    }
    else
    {
        return Vector3D(0, 0, 0);
    }
}