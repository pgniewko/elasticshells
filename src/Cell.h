#ifndef CELL_H
#define	CELL_H

#include <list>
#include <vector>

#include "Environment.h"
#include "force/HookeanForce.h"
#include "force/OsmoticForce.h"
#include "force/HertzianRepulsion.h"
#include "geometry/Vector3D.h"
#include "geometry/Tetrahedron.h"
#include "geometry/Triangle.h"
#include "geometry/Vertex.h"
#include "geometry/VertexTriangle.h"
#include "geometry/algorithms/SimpleTriangulation.h"
#include "simulation/Box.h"
#include "simulation/DomainList.h"

struct cell_params_t
{
    double r_vertex;
    double ecc;
    double dp;
    double gamma;
    double verletR;
    double init_r;
    double vertexVisc;
    double vertexMass;
    double totalVisc;
    double totalMass;
    double growth_rate;
    double vc;
    double bud_d;
    double div_ratio;
};

class Cell
{
    public:
        Cell(int);
        Cell(std::list<Triangle>);
        Cell(const Cell& orig);
        virtual ~Cell();
        double calcSurfaceArea();
        double calcVolume();
        double getMass();
        void calcCM();
        int numberOfTris() ;
        int numberOfVerts();

        void calcBondedForces();
        void calcHarmonicForces();
        void calcOsmoticForces();
        void calcNbForcesON2(const Cell&, Box&);
        void calcNbForcesVL(const Cell&, Box&);
        void calcBoxForces(Box&);

        void voidVerletLsit();
        void builtVerletList(const Cell&, Box&);
        void builtNbList(std::vector<Cell>&, DomainList&, Box&);

        void addVelocity(const Vector3D&);
        void addXYZ(const Vector3D&);

        void setVertexR(double);
        void setEcc(double);
        void setDp(double);
        void setGamma(double);
        void setVisc(double);
        void setMass(double);
        void setCellId(int);
        void setNRT(double);

        void setVerletR(double);
        void setInitR(double);
        
        void setVolumeC(double);
        void setGrowthRate(double);
        void setBudDiameter(double);
        void setDivisionRatio(double);
        
        double getInitR();
        double getVisc();
        Vector3D getCm();
        double getVertexR();
        

        void grow(double);
        void divide();
        
        Vector3D& getVertexXYZ(int);
        Vector3D& getVertexForce(int);

        void voidForces();
        void getDistance(Vector3D&, const Vector3D&, const Vector3D&, Box&);

        Vector3D cm;
        Vertex vertices[MAX_V];
        VertexTriangle triangles[MAX_T];

        int cellId;

    private:
        bool isUnique(std::list<Vector3D>&, const Vector3D&);
        int getVertex(const Vector3D);
        void constructVertices(std::list<Triangle>);
        void constructVTriangles(std::list<Triangle>);
        void constructTopology();
        cell_params_t params;
        int numberV;
        int numberT;
        double nRT;  
};

#endif	/* CELL_H */