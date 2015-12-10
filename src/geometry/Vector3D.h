#ifndef VECTOR3D_H
#define	VECTOR3D_H

#include <iostream>
#include <fstream>
#include <cmath>

#include "Environment.h"
#include "exceptions/DataException.h"

class Vector3D
{
    public:
        inline Vector3D() : x(0), y(0), z(0) {}
        inline Vector3D(double a, double b, double c) : x(a), y(b), z(c) {}
        inline Vector3D(const Vector3D& orig) : x(orig.x), y(orig.y), z(orig.z) {}
        inline virtual ~Vector3D() {};

        inline const Vector3D& operator +=(const Vector3D& v)
        {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }
        
        inline const Vector3D& operator -=(const Vector3D& v)
        {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }
        inline const Vector3D& operator *=(const double a)
        {
            x *= a;
            y *= a;
            z *= a;
            return *this;
        }
        inline const Vector3D& operator /=(const double a)
        {
            return *this *= 1.0 / a;
        }

        inline double length() const
        {
            return sqrt(x * x + y * y + z * z);
        }
        inline double length2() const
        {
            return x * x + y * y + z * z;
        }

        inline void setLength(double r)
        {
            if (length() != 0)
            {
                double rl = r / length();
                x *= rl;
                y *= rl;
                z *= rl;
            }
        }
        
        inline void normalize()
        {
            double len = length();
            x /= len;
            y /= len;
            z /= len;
        }
        
        double x, y, z;
        
        inline double angle(const Vector3D&) const;
};

template <typename InputStreamT>
inline InputStreamT& operator>>(InputStreamT& s, Vector3D& v)
{
    s >> v.x >> v.y >> v.z;
    return s;
}

template <typename OutputStreamT>
inline OutputStreamT& operator<<(OutputStreamT& s, const Vector3D& v)
{
    s <<  v.x << ' ' << v.y << ' ' << v.z << ' ';
    return s;
}

inline Vector3D operator+(const Vector3D& u, const Vector3D& v)
{
    return Vector3D(u.x + v.x, u.y + v.y, u.z + v.z);
}

inline Vector3D operator-(const Vector3D& v)
{
    return Vector3D(-v.x, -v.y, -v.z);
}

inline Vector3D operator -(const Vector3D& u, const Vector3D& v)
{
    return u + -v;
}

inline double operator *(const Vector3D& u, const Vector3D& v)
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

inline bool operator ==(const Vector3D& u, const Vector3D& v)
{
    return (u.x == v.x && u.y == v.y && u.z == v.z);
}

inline bool operator !=(const Vector3D& u, const Vector3D& v)
{
    return (u.x != v.x || u.y != v.y || u.z != v.z);
}

/* Secondary operations */
inline Vector3D operator *(const Vector3D& v, const double a)
{
    return Vector3D(v.x * a, v.y * a, v.z * a);
}

inline Vector3D operator *(const double a, const Vector3D& v)
{
    return v * a;
}

inline Vector3D operator /(const Vector3D& v, const double a)
{
    return v * (1.0 / a);
}

/* Tertiary  operations */
inline Vector3D abs(const Vector3D& v)
{
    return Vector3D(fabs(v.x), fabs(v.y), fabs(v.z));
}

inline Vector3D cross(const Vector3D& v1, const Vector3D& v2)
{
    return Vector3D((v1.y * v2.z - v2.y * v1.z),
                     (v2.x * v1.z - v1.x * v2.z),
                     (v1.x * v2.y - v2.x * v1.y));
}

inline double Vector3D::angle(const Vector3D& v) const
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

//inline double dotd(const Vector3Dd& v1, const Vector3Dd& v2)
//{
//    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
//}

#endif	/* VECTOR3_H */