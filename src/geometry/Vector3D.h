#ifndef VECTOR3D_H
#define	VECTOR3D_H

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

class Vector3D
{
    public:
//        typedef double DataType;

        Vector3D();
        Vector3D(double ix, double iy, double iz);
        Vector3D(const Vector3D& orig);
        virtual ~Vector3D();

        const Vector3D& operator +=(const Vector3D&);
        const Vector3D& operator -=(const Vector3D&);
        const Vector3D& operator *=(const double);
        const Vector3D& operator /=(const double);

        double length() const;
        double length2() const;
        double angle(const Vector3D&) const;
        void setLength(double);

        double x, y, z;

};

template <class InputStreamT>
inline InputStreamT& operator>>(InputStreamT& s, Vector3D& v)
{
    s >> v.x >> v.y >> v.z;
    return s;
}

template <class OutputStreamT>
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

#endif	/* VECTOR3D_H */