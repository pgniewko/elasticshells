#ifndef OSMOTICFORCE_H
#define	OSMOTICFORCE_H

#include "geometry/Vector3D.h"

class OsmoticForce
{
public:
    OsmoticForce();
    OsmoticForce(const OsmoticForce& orig);
    virtual ~OsmoticForce();
    static Vector3D calcForce(const Vector3D& va, const Vector3D& vb, const Vector3D& vc, const Vector3D& vd, const double dp)
    {
        Vector3D BD = vb - vd;
        Vector3D CD = vc - vd;
        Vector3D f = -dp * cross(BD, CD) / 6;
        return f;
    }
private:

};

#endif	/* OSMOTICFORCE_H */

