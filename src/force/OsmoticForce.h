#ifndef OSMOTICFORCE_H
#define	OSMOTICFORCE_H

#include "geometry/Vector3D.h"

class OsmoticForce
{
public:
    OsmoticForce();
    OsmoticForce(const OsmoticForce& orig);
    virtual ~OsmoticForce();
    static Vector3D calcForce(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&, const double);
private:

};

#endif	/* OSMOTICFORCE_H */
