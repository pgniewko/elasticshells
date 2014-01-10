#include "Vector3D.h"

Vector3D::Vector3D() 
{
    x = 0;
    y = 0;
    z = 0;
}

Vector3D::Vector3D(DataType ix, DataType iy, DataType iz)
{
    x = ix; 
    y = iy; 
    z = iz;
}

Vector3D::Vector3D(const Vector3D& orig) 
{
}

Vector3D::~Vector3D() 
{
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
const Vector3D& Vector3D::operator *=(const Vector3D::DataType a)
{
    x *= a;
    y *= a;
    z *= a;
    return *this;
}

const Vector3D& Vector3D::operator /=(const Vector3D::DataType a)
{
    return *this *= 1.0 / a;
}