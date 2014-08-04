#include "Point.h"

//Point::Point() 
//{
//    x = 0.0;
//    y = 0.0;
//    z = 0.0;
//}

Point::Point(float a, float b, float c) : x(a), y(b), z(c) {}

Point::Point(const Point& orig) : x(orig.x), y(orig.y), z(orig.z){}

Point::~Point() {
}

void Point::setLength(float r) 
{
    float rl = r / length();
    x *= rl;
    y *= rl;
    z *= rl;
}

float Point::length() const
{
    return sqrt(x*x + y*y + z*z);
}


Point Point::operator +(const Point& p2) 
{
    return Point(x+p2.x, y+p2.y, z+p2.z);
}

Point Point::operator *(float r) 
{
    return Point(x*r, y*r, z*r);
}

