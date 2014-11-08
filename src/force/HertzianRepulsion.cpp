#include "HertzianRepulsion.h"

HertzianRepulsion::HertzianRepulsion() {}

HertzianRepulsion::HertzianRepulsion(const HertzianRepulsion& orig) {}

HertzianRepulsion::~HertzianRepulsion() {}

Vector3D HertzianRepulsion::calcForce(const Vector3D& dij, const double R, const double E)
{
    return calcForce(dij, R, 0.0, E);
}

Vector3D HertzianRepulsion::calcForce(const Vector3D& dij, const double R1, const double R2, const double E)
{
    double eps;
    double RX;
    
    if (R1 > 0 && R2 > 0)
    {
        RX = R1 * R2 / (R1 + R2);
        eps = R1 + R2 - dij.length();
    }
    else
    {
        RX = std::max(R1, R2);
        eps = RX - dij.length();
    }
    
    if (eps > 0)
    {
        double fmagn = D4_3 * E * pow(RX, 0.5) * pow(eps, 1.5);
        Vector3D f = fmagn * (dij / dij.length()) ;
        return f;
    }
    else
    {
        return Vector3D(0, 0, 0);
    }
}