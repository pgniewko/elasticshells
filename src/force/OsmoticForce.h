#ifndef OSMOTICFORCE_H
#define	OSMOTICFORCE_H

#include "Environment.h"
#include "geometry/Vector3D.h"
#include "geometry/Tetrahedron.h"

class OsmoticForce
{
    public:
        OsmoticForce();
        OsmoticForce(const OsmoticForce& orig);
        virtual ~OsmoticForce();
        static Vector3D calcForce(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&, double, double, const double);
        static void setVolumeFlag(bool);
        static void setEpsilon(double);
        static double getEpsilon();
        static const bool getFlag();


    private:
        static double epsilon;
        static bool volumeFlag;

};

#endif	/* OSMOTICFORCE_H */
