#include "Vector3D.h"

Vector3D::Vector3D() : x(0), y(0), z(0) {}

Vector3D::Vector3D(double a, double b, double c) : x(a), y(b), z(c) {}

Vector3D::Vector3D(const Vector3D& orig) : x(orig.x), y(orig.y), z(orig.z) {}

Vector3D::~Vector3D() {}

double Vector3D::length() const
{
    return sqrt(x * x + y * y + z * z);
}

double Vector3D::length2() const
{
    return x * x + y * y + z * z;
}

void Vector3D::setLength(double r)
{
    if (length() != 0)
    {
        double rl = r / length();
        x *= rl;
        y *= rl;
        z *= rl;
    }
}

double Vector3D::angle(const Vector3D& v) const
{
    if (this->length() == 0 || v.length() == 0)
    {
        throw DataException("Zero length vector");
    }

    double d = *this * v;
    double l1 = length();
    double l2 = v.length();
    double angle = acos(d / (l1 * l2));
    return angle;
}

const Vector3D& Vector3D::operator +=(const Vector3D& v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

const Vector3D& Vector3D::operator -=(const Vector3D& v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

const Vector3D& Vector3D::operator *=(const double a)
{
    x *= a;
    y *= a;
    z *= a;
    return *this;
}

const Vector3D& Vector3D::operator /=(const double a)
{
    return *this *= 1.0 / a;
}