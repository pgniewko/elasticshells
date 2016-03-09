#include "Triangle.h"

//Triangle::Triangle() {}

Triangle::Triangle(Vector3D m, Vector3D n, Vector3D o) : a(m), b(n), c(o) {}


Triangle::Triangle(const Triangle& orig): a(orig.a), b(orig.b), c(orig.c) {}

Triangle::~Triangle() {}

double Triangle::area() const
{
    try
    {
        Vector3D AB = b - a;
        Vector3D AC = c - a;
        double angle = AB.angle(AC);
        double lab = AB.length();
        double lac = AC.length();
        double area = fabs( 0.5 * lab * lac * sin(angle) );
        return area;
    }
    catch (DataException& e)
    {
        return 0.0;
    }
}

double Triangle::min_angle() const
{
    Vector3D ab = a - b;
    Vector3D ac = a - c;
    Vector3D bc = b - c;
    
    double a = ab.length();
    double b = ac.length();
    double c = bc.length();

    double cosA = (a*a + b*b - c*c) / (2*a*b);
    double cosB = (a*a + c*c - b*b) / (2*a*c);
    double cosC = (b*b + c*c - a*a) / (2*b*c);
    
    double thetaA = std::acos(cosA);
    double thetaB = std::acos(cosB);
    double thetaC = std::acos(cosC);
    
    //std::cout << thetaA << " " << thetaB << " " << thetaC << std::endl;
    
    double mintheta = std::numeric_limits<double>::max();
    mintheta = std::min(mintheta, thetaA);
    mintheta = std::min(mintheta, thetaB);
    mintheta = std::min(mintheta, thetaC);
    return mintheta;
}

Vector3D Triangle::normal() const
{
    Vector3D AB = b - a;
    Vector3D AC = c - a;
    Vector3D norm = cross(AB, AC);
    norm /= norm.length();
    return norm;
}

void Triangle::printTriangle() {}

Vector3D Triangle::getVertexA()
{
    return this->a;
}
Vector3D Triangle::getVertexB()
{
    return this->b;
}
Vector3D Triangle::getVertexC()
{
    return this->c;
}