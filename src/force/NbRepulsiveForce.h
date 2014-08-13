#ifndef NBREPULSIVEFORCE_H
#define	NBREPULSIVEFORCE_H

#include "geometry/Vector3D.h"

class NbRepulsiveForce 
{
public:
    NbRepulsiveForce();
    NbRepulsiveForce(const NbRepulsiveForce& orig);
    virtual ~NbRepulsiveForce();
    static Vector3D calcForce(const Vector3D& va, const Vector3D& vb, const double Rc, const double a)
    {
        Vector3D dR = vb - va;
        double r2 = dR.length2();
        
        if (r2 > Rc*Rc) return Vector3D(0,0,0);
       
        double r = dR.length();
        double coeff = ( 1 - Rc / r ) / r2;
        Vector3D f  = coeff * dR;
        return f;
    }
private:

};

#endif	/* NBREPULSIVEFORCE_H */
