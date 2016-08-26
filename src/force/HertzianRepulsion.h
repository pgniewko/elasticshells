#ifndef HERTZIANREPULSION_H
#define	HERTZIANREPULSION_H

#include <cmath>

#include "Environment.h"
#include "geometry/Vector3D.h"

class HertzianRepulsion
{
    public:
        HertzianRepulsion() = delete;
        HertzianRepulsion(const HertzianRepulsion& orig) = delete;
        virtual ~HertzianRepulsion() = delete;

        static Vector3D calcForce(const Vector3D&, const double, const double, const double, const double, const double, const double);
        static double calcEnergy(const Vector3D&, const double, const double, const double, const double, const double, const double);
    private:

};

#endif	/* HERTZIANREPULSION_H */

