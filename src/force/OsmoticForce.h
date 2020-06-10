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
        static void set_volume_flag(bool);
        static void set_epsilon(double);
        static double get_epsilon();
        static bool get_flag();

    private:
        static double epsilon;
        static bool volume_flag;
};

#endif	/* OSMOTICFORCE_H */
