#ifndef HERTZIANREPULSION_H
//#ifndef BOXCELLREPULSION_H
#define	HERTZIANREPULSION_H

#include <math.h>

#include "Environment.h"
#include "geometry/Vector3D.h"

class HertzianRepulsion {
public:
    HertzianRepulsion();
    HertzianRepulsion(const HertzianRepulsion& orig);
    virtual ~HertzianRepulsion();
    
    static Vector3D calcForce(Vector3D&, double, double);
private:

};

#endif	/* HERTZIANREPULSION_H */

