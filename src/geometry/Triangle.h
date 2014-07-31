#ifndef TRIANGLE_H
#define	TRIANGLE_H

#include "Point.h"

class Triangle {
public:
//    Triangle();
    Triangle(Point, Point, Point);
    Triangle(const Triangle& orig);
    virtual ~Triangle();
    Point a, b, c;
     
private:
   

};

#endif	/* TRIANGLE_H */

