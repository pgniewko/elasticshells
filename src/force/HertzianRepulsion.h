#ifndef HERTZIANREPULSION_H
#define	HERTZIANREPULSION_H

#include <cmath>

#include "Environment.h"
#include "geometry/Vector3D.h"

class HertzianRepulsion
{
    public:
        HertzianRepulsion();
        HertzianRepulsion(const HertzianRepulsion& orig);
        virtual ~HertzianRepulsion();

        //static Vector3D calcForce(const Vector3D&, const double, const double);
        static Vector3D calcForce(const Vector3D&, const double, const double, const double, const double, const double, const double);
    private:

};

#endif	/* HERTZIANREPULSION_H */

