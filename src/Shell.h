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

struct shell_params_t
{
    double vertex_r;
    double ecc;
    double nu;
    double dp;
    double init_r;
    double vol_c;
};


class Shell
{
        friend class Tinker;
        friend class DomainList;
        friend class Energy;
        friend class Restarter;

    public:

        Shell();
        Shell(std::list<Triangle>);
        Shell(const Shell& orig);
        virtual ~Shell();
        double calcSurfaceArea() const;
        double calcVolume(double = 0.0) const;
        void calc_cm();
        int getNumberTriangles() const;
        int getNumberVertices() const;
        int getNumberHinges() const;
        
        void add_vector(const Vector3D&);

        void setVertexR(double);
        void setEcc(double);
        void setDp(double);
        void pushDp(double dp_)
        {
            params.dp = dp_;
        }
        void setDp(double, double);
        void setSpringConst(double, double, double);
        void setShellId(int);
        void setNu(double);
        void setBSprings(double, double, double);
        void setConstantVolume(double = 1.0);
        double checkVolumeCondition();
        double ajust_turgor(double = 0.0);

        void setInitR(double);

        double getInitR() const;
        Vector3D getCm() const;
        double getVertexR() const;
        double getE() const;
        double getNu() const;
        
        Vector3D center_of_mass;
        
        std::vector<Vertex> vertices;
        std::vector<VertexTriangle> triangles;
        std::vector<BendingHinge> bhinges;
        
        int shell_id = -1;
        double calcSurfaceArea(double) const;

        double getTurgor() const;
        void update(double = 0.0);

        const shell_params_t& get_params() const;

        static bool bending;

        friend std::ostream& operator<< (std::ostream&, const Shell&);

        static double FORCE_FRAC;
        static double MIN_FORCE;

    private:

        void randomRotate();
        
        shell_params_t params;
        int number_v = 0;
        int number_t = 0;
        int number_s = 0;
        double nRT = 0.0;
        double V0 = 0.0;

        static utils::Logger cell_log;
};

#endif	/* CELL_H */
