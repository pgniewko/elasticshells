#include "NbRepulsiveForce.h"

NbRepulsiveForce::NbRepulsiveForce() {}

NbRepulsiveForce::NbRepulsiveForce(const NbRepulsiveForce& orig) {}

NbRepulsiveForce::~NbRepulsiveForce() {}

/* dR must point to the vertex to which force is added
 */
Vector3D NbRepulsiveForce::calcForce(const Vector3D& dR, const double Rc, const double a)
{
    double r2 = dR.length2();

    if (r2 > Rc * Rc)
    {
        return Vector3D(0, 0, 0);
    }

    double sigr2 = Rc * Rc / r2;
    double r6 = sigr2 * sigr2 * sigr2;
    double r12 = r6 * r6;
    Vector3D f = 12 * a * (r12 - r6) * dR / r2; // - 12 * a * (r^-12 - r^-6) / r * versor{R};  versor{R} = dR / |dR|
    return f;   
}

Vector3D NbRepulsiveForce::calcForce(const Vector3D& va, const Vector3D& vb, const double Rc, const double a)
{  
    Vector3D dR = va - vb;
    Vector3D res = calcForce(dR, Rc, a);
    return res;
}