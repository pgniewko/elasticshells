#ifndef HOOKEANFORCE_H
#define	HOOKEANFORCE_H

#include "geometry/Vector3D.h"

class HookeanForce
{
public:
    HookeanForce();
    HookeanForce(const HookeanForce& orig);
    virtual ~HookeanForce();
    
    static Vector3D calcForce(const Vector3D&, const Vector3D&, const double, const double);

};
#endif	/* HOOKEANFORCE_H */
