#ifndef VERTEXTRIANGLE_H
#define	VERTEXTRIANGLE_H

#include "Environment.h"
#include "Vertex.h"
#include "Triangle.h"
#include "geometry/Vector3D.h"

class Vertex;

class VertexTriangle
{
    public:
        VertexTriangle();
        VertexTriangle(int, int, int);
        VertexTriangle(const VertexTriangle& orig);
        virtual ~VertexTriangle();
        void setId(int);
        int getId() const;
        double area(const Vertex[]) const;
        double area(const Vertex[], const Vector3D, double) const;

        Vector3D normal(const Vertex[]) const;

        void printVertexTriangle() const;
        void subsVertex(int, int);

        int ia = -1;
        int ib = -1;
        int ic = -1;
        int myindex = -1;
};

#endif	/* VERTEXTRIANGLE_H */
