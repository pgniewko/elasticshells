#ifndef BENDINGSPRING_H
#define	BENDINGSPRING_H

#include <iostream>
#include <iomanip>

#include "Environment.h"
#include "Vector3D.h"
#include "Triangle.h"
#include "Vertex.h"

class BendingSpring {
public:
    BendingSpring();
    BendingSpring(int, int, int, int);
    BendingSpring(const BendingSpring& orig);
    virtual ~BendingSpring();
    
    void setD(const double&, const double&, const double&);
    void calcBendingForces(Vertex[]) const;
    void setThetaZero(const Vertex[]);
    friend bool operator== (BendingSpring&, BendingSpring&);
    
private:
    double D = 0.0;
    double sinHalfTheta0 = 0.0;
    int x1 = -1;
    int x2 = -1;
    int x3 = -1;
    int x4 = -1;
    double calcSinTheta(const Vertex[]) const;
};

inline bool operator== (BendingSpring& bs1, BendingSpring& bs2)
{
    return (bs1.x1 == bs2.x1 && bs1.x2 == bs2.x2 && bs1.x3 == bs2.x3 && bs1.x4 == bs2.x4);
}

#endif	/* BENDINGSPRING_H */

