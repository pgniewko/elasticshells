#ifndef HOOKEANFORCE_H
#define	HOOKEANFORCE_H

#include "geometry/Vector3D.h"

class HookeanForce
{
public:
    HookeanForce();
    HookeanForce(const HookeanForce& orig);
    virtual ~HookeanForce();
    
    static Vector3D calcForce(const Vector3D& va, const Vector3D& vb, const double R0, const double gamma)
    {
        Vector3D dR = vb - va;
        double R = dR.length();
        Vector3D f = gamma * dR * (1 - R0 / R);
        return f;
    }

};
#endif	/* HOOKEANFORCE_H */

