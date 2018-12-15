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
#include "geometry/Element.h"
#include "geometry/Hinge.h"
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
        friend class Energy;
        friend class Restarter;

    public:

        Shell();
        Shell(int,int,int);
        Shell(std::list<Triangle>);
        Shell(const Shell& orig);
        virtual ~Shell();
        double calc_volume(double = 0.0) const;
        double calcSurfaceArea() const;
        double calcSurfaceArea(double) const;
        void calc_cm();
        int getNumberTriangles() const;
        int getNumberVertices() const;
        int getNumberHinges() const;
        
        void add_vector(const Vector3D&);

        void set_vertex_size(double);
        void set_ecc(double);
        void set_dp(double);
        void push_dp(double dp_)
        {
            params.dp = dp_;
        }
        void set_dp(double, double);
        void set_elements_parameters(double, double, double);
        void set_shell_id(int);
        void set_nu(double);
        void set_hinges(double, double, double);
        void set_constant_volume(double = 1.0);
        double check_volume_condition();
        double ajust_turgor(double = 0.0);

        void setInitR(double);

        double get_r0() const;
        Vector3D get_cm() const;
        double get_vertex_size() const;
        double get_E() const;
        double get_nu() const;
        
        Vector3D center_of_mass;
        
        std::vector<Vertex> vertices;
        std::vector<Element> triangles;
        std::vector<Hinge> hinges;
        
        int shell_id = -1;

        double get_turgor() const;
//        void update();

        const shell_params_t& get_params() const;

        static bool bending;

        friend std::ostream& operator<< (std::ostream&, const Shell&);

        static double FORCE_FRAC;
        static double MIN_FORCE;

    private:

        void random_rotate();
        
        shell_params_t params;
        int number_v = 0;
        int number_t = 0;
        int number_h = 0;
        double nRT = 0.0;
        double V0 = 0.0;

        static utils::Logger cell_log;
};

#endif	/* CELL_H */
