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
        Vector3D normal() const;
        void printTriangle();
        Vector3D a, b, c;

        Vector3D getVertexA();
        Vector3D getVertexB();
        Vector3D getVertexC();

};

#endif	/* TRIANGLE_H */
