#include "Triangle.h"

//Triangle::Triangle() 
//{
//}

//Triangle::Triangle(Vector3D m, Vector3D n, Vector3D o) : a(m), b(n), c(o) {}

Triangle::Triangle(Vector3D& m, Vector3D& n, Vector3D& o)
{
    a = m;
    b = n;
    c = o;
    cout << "Address of m : "<< m << endl;
    cout << "Address of a : "<< a << endl;
}

Triangle::Triangle(const Triangle& orig): a(orig.a), b(orig.b), c(orig.c) {}

Triangle::~Triangle() {}

double Triangle::area() const
{
    return 0.0;
//    Vector3D AB = b - a;
//    Vector3D AC = c - a;
//    double angle = AB.angle(AC); 
//    double lab = AB.length();
//    double lac = AC.length();
//    double area = fabs( 0.5 * lab * lac * sin(angle) );
//    return area;
}

void Triangle::setId(int idx)
{
    id = idx;
}

int Triangle::getId()
{
    return id;
}

void Triangle::printTriangle()
{
    cout << "mem of a" << &a;
    cout << " " << a;
    cout << "mem of b" << &b;
    cout << " " << b;
    cout << "mem of c" << &c;
    cout << " " << c;
    
    cout << endl ;
}
