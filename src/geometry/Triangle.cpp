#include "Triangle.h"

//Triangle::Triangle() 
//{
//}

Triangle::Triangle(Point m, Point n, Point o) : a(m), b(n), c(o) { }

Triangle::Triangle(const Triangle& orig): a(orig.a), b(orig.b), c(orig.c) {}

Triangle::~Triangle() {
}

