#include "Triangle.h"

//Triangle::Triangle() 
//{
//}

Triangle::Triangle(Vector3D m, Vector3D n, Vector3D o) : a(m), b(n), c(o) {}

Triangle::Triangle(const Triangle& orig): a(orig.a), b(orig.b), c(orig.c) {}

Triangle::~Triangle() {}

