#ifndef TRIANGLE_H
#define	TRIANGLE_H

#include "Environment.h"
#include "Vector3D.h"

class Triangle
{
    public:
        Triangle(Vector3D, Vector3D, Vector3D);
        Triangle(const Triangle& orig);
        virtual ~Triangle();
        double area() const;
        double min_angle() const;
        Vector3D normal() const;
        Vector3D a, b, c;
};

#endif	/* TRIANGLE_H */
