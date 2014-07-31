#ifndef POINT_H
#define	POINT_H

#include <fstream>
#include <cmath>

class Point {
public:
//    Point();
    Point(float, float, float);
    Point(const Point& orig);
    virtual ~Point();
    
    Point operator +(const Point& p);
    Point operator *(float r);
    float length() const;
    void setLength(float r);
    
    float x, y, z;
    
private:

};

#endif	/* POINT_H */

