#include "HertzianRepulsion.h"

Vector3D HertzianRepulsion::calcForce(const Vector3D& dij, const double R1, const double R2, const double e1, const double e2, const double nu1, const double nu2)
{
    double h;
    double r_eff;
    double e_eff;

    if (R1 > 0 && R2 > 0)
    {
        r_eff = R1 * R2 / (R1 + R2);
        h = R1 + R2 - dij.length();
    }
    else
    {
        r_eff = std::max(R1, R2);
        h = r_eff - dij.length();
    }

    e_eff = (1 - nu1 * nu1) / e1 + (1 - nu2 * nu2) / e2;
    e_eff = 1.0 / e_eff;

    if (h > 0)
    {
        double fmagn = constants::d4_3 * e_eff * pow(r_eff, 0.5) * pow(h, 1.5);
        Vector3D f = fmagn * (dij / dij.length()) ;
        return f;
    }
    else
    {
        return Vector3D(0, 0, 0);
    }
}

double HertzianRepulsion::calcEnergy(const Vector3D& dij, const double R1, const double R2, const double e1, const double e2, const double nu1, const double nu2)
{
    double energy = 0.0;
    double h;
    double r_eff;

    if (R1 > 0 && R2 > 0)
    {
        r_eff = R1 * R2 / (R1 + R2);
        h = R1 + R2 - dij.length();
    }
    else
    {
        r_eff = std::max(R1, R2);
        h = r_eff - dij.length();
    }



    double D = 0.75 * ( (1 - nu1 * nu1) / e1 + (1 - nu2 * nu2) / e2 );

    if (h > 0)
    {
        energy = constants::d2_5 * pow(r_eff, 0.5) * pow(h, 2.5) / D;
        return energy;
    }
    else
    {
        return 0.0;
    }
}