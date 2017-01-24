#ifndef PACKER_H
#define PACKER_H

// TODO: at some point this class should become an independent library.

#include "Environment.h"
#include "geometry/Vector3D.h"
#include "Cell.h"
#include "simulation/Box.h"
#include "force/HertzianRepulsion.h"
#include "utils/io/LogSimulation.h"

struct point_t
{
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

class Packer {
public:
    Packer();
    Packer(const Packer& orig);
    virtual ~Packer();
    
    static void packCells(Box&, std::vector<Cell>&, double);
    
private:
    
    static void fire(std::vector<point_t>&, box_t&);
    static void velocityVerlet(std::vector<point_t>&, box_t&);
    static void calcForces(std::vector<point_t>&, box_t&);
    static bool check_min_force(std::vector<point_t>&, box_t&);
    
    static bool jammed(std::vector<point_t>&, box_t&);
    static void inflatePoints(std::vector<point_t>&);
    static void recenterCells(std::vector<point_t>&, box_t&);
    
    static double calcPressure(std::vector<point_t>&, box_t&);
    static double boxForce(point_t&, box_t&);
    
    static void getDistance(Vector3D&, const Vector3D&, const Vector3D&, const box_t&);
    static void calcBoxForces(std::vector<point_t>&, const box_t&);
    
    
    static bool anyRattlers(std::vector<point_t>&, box_t&);
    static int boxContacts(point_t&, box_t&);
    static int cellContacts(point_t&, point_t&, box_t&);
    
    
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

