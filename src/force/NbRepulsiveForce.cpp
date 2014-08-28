#include "NbRepulsiveForce.h"

NbRepulsiveForce::NbRepulsiveForce() {}

NbRepulsiveForce::NbRepulsiveForce(const NbRepulsiveForce& orig) {}

NbRepulsiveForce::~NbRepulsiveForce() {}

static Vector3D NbRepulsiveForce::calcForce(const Vector3D& va, const Vector3D& vb, const double Rc, const double a)
{
        Vector3D dR = vb - va;
        double r2 = dR.length2();
        
        if (r2 > Rc*Rc) return Vector3D(0, 0, 0);
       
        double r = dR.length();
        double coeff = a * ( 1 - Rc / r ) / r2;
        Vector3D f = coeff * dR;
        return f;
}