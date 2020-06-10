#ifndef PACKER_H
#define PACKER_H

// TODO: at some point this class should become an independent library.

#include "Environment.h"
#include "geometry/Vector3D.h"
#include "Shell.h"
#include "simulation/Box.h"
#include "force/HertzianRepulsion.h"
#include "utils/io/LogSimulation.h"

struct point_t
{
    point_t() :radius(0), radius_f(0)  {}
    Vector3D r_c;
    Vector3D v_c;
    Vector3D f_c;
    Vector3D f_p;
    double radius;
    double radius_f;
};

struct box_t
{
    double x;
    double y;
    double z;
    bool pbc;
};

class Packer
{
    public:
        Packer();
        Packer(const Packer& orig);
        virtual ~Packer();

        static void pack_shells(Box&, std::vector<Shell>&, double, bool = true);

    private:

        static void fire(std::vector<point_t>&, box_t&);
        static void velocity_verlet(std::vector<point_t>&, box_t&);
        static void calc_forces(std::vector<point_t>&, box_t&);
        static bool check_min_force(std::vector<point_t>&, box_t&);

        static bool jammed(std::vector<point_t>&, box_t&, double&);
        static void inflate_points(std::vector<point_t>&);
        static void recenter_shells(std::vector<point_t>&, box_t&);

        static double calc_pressure(std::vector<point_t>&, box_t&);
        static double box_force(point_t&, box_t&);

        static void get_distance(Vector3D&, const Vector3D&, const Vector3D&, const box_t&);
        static void calc_box_forces(std::vector<point_t>&, const box_t&);


        static bool any_rattler_or_crowder(const std::vector<point_t>&, const box_t&, double&);
        static int box_contacts(const point_t&, const box_t&);
        static int shell_contacts(const point_t&, const point_t&, const box_t&);

        static int FIRE_Nmin;
        static int FIRE_N;
        static double FIRE_DT;
        static double FIRE_ALPHA;
        static double FIRE_DTMAX;

        static double MIN_FORCE;
        static double r_ext;

        static double P_MIN;
        static double P_MAX;

        static utils::Logger packer_logs;

};

#endif /* PACKER_H */

