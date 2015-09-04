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
        int getId();
        double area(const Vertex[]);
        double area(const Vertex[], const Vector3D, double);

        Vector3D normal(const Vertex[]);

        void printVertexTriangle();
        void subsVertex(int, int);

        int ia, ib, ic;
        int myindex;
};

#endif	/* VERTEXTRIANGLE_H */
