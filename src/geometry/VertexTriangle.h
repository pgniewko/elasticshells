#ifndef VERTEXTRIANGLE_H
#define	VERTEXTRIANGLE_H

#include "Environment.h"
#include "geometry/Triangle.h"
#include "geometry/Vector3D.h"
#include "geometry/Vertex.h"

class Vertex;

class VertexTriangle
{
        friend class Restarter;
    public:
        VertexTriangle();
        explicit VertexTriangle(int, int, int);
        VertexTriangle(const VertexTriangle& orig);
        virtual ~VertexTriangle();
        void setId(int);
        int getId() const;
        double area(const Vertex[]) const;
        double area(const Vertex[], const Vector3D, double) const;

        Vector3D normal(const Vertex[]) const;
        void calcFemForces(Vertex[]) const;
        double calcFemEnergy(const Vertex[]) const;

        void printVertexTriangle() const;
        void subsVertex(int, int);

        void setParams(const Vertex[], const double, const double, const double);

        int ia = -1;
        int ib = -1;
        int ic = -1;
        int myid = -1;

        friend std::ostream& operator<< (std::ostream&, const VertexTriangle&);

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
