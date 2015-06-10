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
#include "simulation/Tinker.h"

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

enum class cell_phase_t
{
    C_G0,   // gap phase (resting)
    C_G1,   // mother cell growth
    C_SG2,  // S+G2 - i.e. bud creation and budding phase
    C_M     // cell division phase
};

class Cell
{
        friend class Tinker;
    public:

        Cell(int);
        Cell(std::list<Triangle>);
        Cell(const Cell& orig);
        virtual ~Cell();
        double calcSurfaceArea();
        double calcVolume();
        double calcVolume(double);
        double getMass();
        void calcCM();
        int numberOfTris() ;
        int numberOfVerts();

        //void calcR0av();

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
        //void moveToXYZ(const Vector3D&);

        void setVertexR(double);
        void setEcc(double);
        void setDp(double);
        void setGamma(double);
        void setVisc(double);
        void setMass(double);
        void setCellId(int);
        void setNRT(double);
        void setNRT(double, double);

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

        Vector3D& getVertexXYZ(int);
        Vector3D& getVertexForce(int);

        void voidForces();
        void getDistance(Vector3D&, const Vector3D&, const Vector3D&, Box&);

        void cellCycle();

        Vector3D cm_m;
        Vector3D cm_b;
        Vertex vertices[MAX_V];
        VertexTriangle triangles[MAX_T];

        int cellId;

        void setR0AvForAll();
        double getPercLength(int, int);
        double nbMagnitudeForce(std::vector<Cell>, Box&, int);
        double nbMagnitudeForce(Cell, Box&);

        double contactForce(const Cell&, Box&);
        //double contactForce(Box&);

        double contactArea(const Cell&, Box&);
        double contactArea(Box&);

        double surfaceStrainEnergy();

    private:
        void grow();
        void bud();
        void divide();
        void findBud();
        void randomRotate();

        void calcAverageR0();
        double sumL2();
        cell_params_t params;
        cell_phase_t my_phase;
        int numberV;
        int numberT;
        double nRT;

        double r0av;
        int vcounter; // remove - prawdopodobnie bezurzyteczny

        int budIdx[MAX_V];
        int budVno;

        static utils::Logger cell_log;
};


#endif	/* CELL_H */