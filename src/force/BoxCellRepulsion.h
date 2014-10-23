#ifndef BOXCELLREPULSION_H
#define	BOXCELLREPULSION_H

#include <math.h>

#include "Environment.h"
#include "geometry/Vector3D.h"

class BoxCellRepulsion {
public:
    BoxCellRepulsion();
    BoxCellRepulsion(const BoxCellRepulsion& orig);
    virtual ~BoxCellRepulsion();
    
    static Vector3D calcForce(Vector3D&, double, double);
private:

};

#endif	/* BOXCELLREPULSION_H */

