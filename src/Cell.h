#ifndef CELL_H
#define	CELL_H

#include <list>
#include <vector>

#include "Environment.h"
#include "force/HookeanForce.h"
#include "force/OsmoticForce.h"
#include "force/NbRepulsiveForce.h"
#include "force/BoxCellRepulsion.h"
#include "geometry/Vector3D.h"
#include "geometry/Tetrahedron.h"
#include "geometry/Triangle.h"
#include "geometry/Vertex.h"
#include "geometry/VertexTriangle.h"
#include "geometry/algorithms/SimpleTriangulation.h"
#include "simulation/Box.h"
#include "simulation/DomainList.h"

class Cell
{
    public:
        Cell(int);
        Cell(list<Triangle>);
        Cell(const Cell& orig);
        virtual ~Cell();
        double calcSurfaceArea();
        double calcVolume();
        double getMass();
        void calcCM();
        int numberofFaces() ;
        int numberofVertices();
        void builtVerletList(const Cell&, Box&);
        void voidVerletLsit();
        
        void builtNbList(vector<Cell>&, DomainList&, Box&);

        void calcForces();
        void calcForces(const Cell&, Box&);
        void calcForcesVL(const Cell&, Box&);
        void calcForces(Box&);
        //double calcBoxForces(Box&); // TODO: refactor it!
        //void calcStressTensor(Box&, double*); //TODO: refacoting !
        
        void addVelocity(const Vector3D&);
        void addXYZ(const Vector3D&);

        void setRc(double);
        void setA(double);
        void setDp(double);
        void setGamma(double);
        void setVisc(double);
        void setMass(double);
        void setCellId(int);
        void setRCellBox(double);
        void setNRT(double);

        void setVerletR(double);
        void setInitR(double);

        double getInitR();
        double getVisc();
        Vector3D getCm();
        double getRcb();
        
        Vector3D getVertexXYZ(int);
        Vector3D getVertexForce(int);

        void voidForces();
        void getDistance(Vector3D& djk, const Vector3D& vj, const Vector3D& vk, Box&);

        Vector3D cm;
        Vertex vertices[MAX_V];
        VertexTriangle triangles[MAX_T];

        int cellId;

    private:
        bool isUnique(list<Vector3D>&, const Vector3D&);
        int getVertex(const Vector3D);
        void constructVertices(list<Triangle>);
        void constructVTriangles(list<Triangle>);
        void constructTopology();
        int numberV;
        int numberT;
        double Rc;
        double rCellBox;
        double a;
        double dp;
        double gamma;
        double verletR;
        double initR;
        double visc0;
        double mass0;
        double visc0tot;
        double mass0tot;
        
        double nRT;
};

#endif	/* CELL_H */