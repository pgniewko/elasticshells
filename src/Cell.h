#ifndef CELL_H
#define	CELL_H

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
#include "force/BendingHinge.h"
#include "geometry/algorithms/SimpleTriangulation.h"
#include "simulation/Box.h"
#include "simulation/DomainList.h"
#include "simulation/Tinker.h"

struct cell_params_t
{
    double vertex_r;
    double ecc;
    double nu;
    double dp;
    double verlet_f;
    double init_r;
};

class Cell
{
        friend class Tinker;
    public:

        Cell(int);
        Cell(std::list<Triangle>);
        Cell(const Cell& orig);
        virtual ~Cell();
        double calcSurfaceArea() const;
        double calcVolume(double = 0.0) const;
        void calcCM();
        int getNumberTriangles() const;
        int getNumberVertices() const;

        void calcBondedForces();
        void calcHarmonicForces();
        void calcFemForces();
        void calcOsmoticForces();
        void calcNbForcesON2(const Cell&, const Box&);
        void calcNbForcesVL(const Cell&, const Box&);
        void calcBoxForces(const Box&);

        void voidVerletLsit();
        void builtVerletList(const Cell&, const Box&);
        void builtNbList(std::vector<Cell>&, DomainList&, const Box&);

        void addXYZ(const Vector3D&);

        void setVertexR(double);
        void setEcc(double);
        void setDp(double);
        void setDp(double, double);
        void setSpringConst(double, double, double, char*);
        void setCellId(int);
        void setNu(double);
        void setBSprings(double, double, double);

        void setVerletR(double);
        void setInitR(double);

        double getInitR() const;
        Vector3D getCm() const;
        double getVertexR() const;
        double getE() const;
        double getNu() const;

        void voidForces();
        void getDistance(Vector3D&, const Vector3D&, const Vector3D&, const Box&) const;

        Vector3D cm_m;
        Vertex vertices[MAX_V];
        VertexTriangle triangles[MAX_T];
        BendingHinge bhinges [2 * MAX_T];

        int cell_id = -1;
        double contactForce(const Cell&, const Box&) const;
        double contactForce(const Box&) const;
        double contactForceSF(const Box&) const; // for Surface Force use
        double contactArea(const Cell&, const Box&) const;
        double contactArea(const Box&, double = 0.0) const;
        double activeArea(const Box&, const std::vector<Cell>&, double&, bool = false) const;
        double activeAreaFraction(const Box&, const std::vector<Cell>&, double&, bool = false) const;
        double strainEnergy(const Box&) const;
        double maxStrain() const;
        double minStrain() const;
        double nbIntra(const Box&) const;
        double getTurgor() const;
        double getStrain(int, int) const;
        void update(double = 0.0);

        static bool no_bending;

    private:

        void randomRotate();

        bool isInContact(const int, const Cell&, const Box&) const;
        bool isInContact(const int, const Box&) const;
        double project_force(const Cell&, const Box&, const Vector3D&, const int) const;
        double project_force(const Box&, const Vector3D&, const int) const;
        Vector3D box_force(const Box&, const int) const;

        double sumL2() const;

        cell_params_t params;
        int number_v = 0;
        int number_t = 0;
        int number_s = 0;
        double nRT = 0.0;
        double V0 = 0.0;

        bool fem_flag = false;
        bool bending_flag = true;

        static utils::Logger cell_log;
};


#endif	/* CELL_H */
