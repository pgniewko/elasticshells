#include "NbRepulsiveForce.h"

NbRepulsiveForce::NbRepulsiveForce() {}

NbRepulsiveForce::NbRepulsiveForce(const NbRepulsiveForce& orig) {}

NbRepulsiveForce::~NbRepulsiveForce() {}

Vector3D NbRepulsiveForce::calcForce(const Vector3D& va, const Vector3D& vb, const double Rc, const double a)
{
        Vector3D dR = vb - va;
        double r2 = dR.length2();
        
        if (r2 > Rc * Rc) return Vector3D(0, 0, 0);
       
        //double r = Rc / dR.length();
        double sigr2 = Rc * Rc / r2;
        double r6 = sigr2 * sigr2 * sigr2;
        double r12 = r6 * r6;
        double coeff = -a / r2;
        Vector3D f = coeff * dR * (r12 - r6); // - 12 * a * (r^-12 - r^-6) / r * versor{R};  versor{R} = dR / |dR|
        return f;
        
        //double r = dR.length();
        //double coeff = a * ( 1 - Rc / r ) / r2;
        //Vector3D f = coeff * dR;
        //return f;
}