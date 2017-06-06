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
        
        static Vector3D calcForce(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&, const double);
        static void setVolumeFlag(bool);
        static void setEpsilon(double);
        static double getEpsilon();
        static bool getFlag();


    private:
        static double epsilon;
        static bool volumeFlag;
};


inline Vector3D OsmoticForce::calcForce(const Vector3D& va, const Vector3D& vb, const Vector3D& vc, const Vector3D& vd, const double turgor)
{
    Vector3D BD = vb - vd;
    Vector3D CD = vc - vd;
    Vector3D f = cross(BD, CD) / 6;
    f *= turgor;

    return f * Tetrahedron::volumeSgn(va, vb, vc, vd);
}

#endif	/* OSMOTICFORCE_H */
