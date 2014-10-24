#include "BoxCellRepulsion.h"

BoxCellRepulsion::BoxCellRepulsion() {}

BoxCellRepulsion::BoxCellRepulsion(const BoxCellRepulsion& orig) {}

BoxCellRepulsion::~BoxCellRepulsion() {}

Vector3D BoxCellRepulsion::calcForce(Vector3D& dij, double E, double R)
{
    double eps = R - dij.length();
    if (eps > 0)
    {
        double fmagn = D4_3 * E * pow(R, 0.5) * pow(eps, 1.5) / dij.length();
        Vector3D f = fmagn * dij ;
        return f;
    }
    else
    {
        return Vector3D(0, 0, 0);
    }
}