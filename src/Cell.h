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
#include "simulation/Tinker.h"
//#include "simulation/Restarter.h"

struct cell_params_t
{
    double vertex_r;
    double ecc;
    double nu;
    double dp;
    double init_r;
    double vol_c;
};

class Cell
{
        friend class Tinker;
        friend class DomainList;
        friend class Energy;
        friend class Restarter;
    public:

        Cell();
        explicit Cell(int);
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
        void calcBoxForces(const Box&);

        void addXYZ(const Vector3D&);

        void setVertexR(double);
        void setEcc(double);
        void setDp(double);
        void pushDp(double dp_) {params.dp = dp_;}
        void setDp(double, double);
        void setSpringConst(double, double, double, std::string);
        void setCellId(int);
        void setNu(double);
        void setBSprings(double, double, double);
        void setConstantVolume(double = 1.0);
        double checkVolumeCondition(double = 0.0);
        void ajustTurgor(double = 0.0);

        void setInitR(double);

        double getInitR() const;
        Vector3D getCm() const;
        double getVertexR() const;
        double getE() const;
        double getNu() const;

        void voidForces();

        Vector3D cm_m;
        Vertex vertices[MAX_V];
        VertexTriangle triangles[MAX_T];
        BendingHinge bhinges [2 * MAX_T];

        int cell_id = -1;
        double contactForce(const Cell&, const Box&) const;
        double contactForce(const Box&) const;
        double contactForceSF(const Box&) const; // for Surface Force use
        
        double contactArea(const Cell&, const Box&) const;
        double contactArea(const Cell&, const Box&, const double) const;
        double contactArea(const Box&, double = 0.0) const;
        
        
        double activeArea(const Box&, const std::vector<Cell>&, double = 0.0) const;
        double calcSurfaceArea(double) const;
        double contactArea2(const Box&, double = 0.0) const;
        //double activeAreaFraction(const Box&, const std::vector<Cell>&, double&, double = 0.0) const;
        
        double strainEnergy(const Box&) const;
        double maxStrain() const;
        double minStrain() const;
        double nbIntra(const Box&) const;
        double getTurgor() const;
        double getStrain(int, int) const;
        void update(double = 0.0);

        const cell_params_t& get_params() const;

        static bool no_bending;

        friend std::ostream& operator<< (std::ostream&, const Cell&);

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
        //bool bending_flag = true;

        static utils::Logger cell_log;
};

#endif	/* CELL_H */
