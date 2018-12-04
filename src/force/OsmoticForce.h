#ifndef OSMOTICFORCE_H
#define	OSMOTICFORCE_H

#include "Environment.h"
#include "geometry/Vector3D.h"
#include "geometry/Tetrahedron.h"

class OsmoticForce
{
    public:
        OsmoticForce() = delete;
        OsmoticForce(const OsmoticForce& orig) = delete;
        virtual ~OsmoticForce() = delete;
        static void setVolumeFlag(bool);
        static void setEpsilon(double);
        static double getEpsilon();
        static bool getFlag();


    private:
        static double epsilon;
        static bool volumeFlag;
};

#endif	/* OSMOTICFORCE_H */
