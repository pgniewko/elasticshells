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
        void calcFemForces(Vertex[]) const;

        void printVertexTriangle() const;
        void subsVertex(int, int);

        void setParams(const Vertex[], const double, const double, const double);

        int ia = -1;
        int ib = -1;
        int ic = -1;
        int myindex = -1;

    private:
        void setL2(const Vertex[]);
        void setAn(const Vertex[]);
        void setKi(const Vertex[], const double&, const double&, const double&);
        void setCi(const Vertex[], const double&, const double&, const double&);

        double an[3];
        double L2[3];
        double ki[3];
        double ci[3];

};

#endif	/* VERTEXTRIANGLE_H */
