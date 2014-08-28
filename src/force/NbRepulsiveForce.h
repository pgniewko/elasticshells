#ifndef NBREPULSIVEFORCE_H
#define	NBREPULSIVEFORCE_H

#include "geometry/Vector3D.h"

class NbRepulsiveForce 
{
public:
    NbRepulsiveForce();
    NbRepulsiveForce(const NbRepulsiveForce& orig);
    virtual ~NbRepulsiveForce();
    static Vector3D calcForce(const Vector3D&, const Vector3D&, const double, const double);

private:

};

#endif	/* NBREPULSIVEFORCE_H */
