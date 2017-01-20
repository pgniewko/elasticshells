#ifndef PACKER_H
#define PACKER_H

#include "Environment.h"
#include "geometry/Vector3D.h"
#include "Cell.h"
#include "simulation/Box.h"
#include "force/HertzianRepulsion.h"

struct point_t
{
    Vector3D r_c;
    Vector3D v_c;
    Vector3D f_c;
    Vector3D f_p;
    double radius;
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
    
    static int FIRE_Nmin;
    static int FIRE_N;
    static double FIRE_DT;
    static double FIRE_ALPHA;
    static double FIRE_DTMAX;
    
    static double MIN_FORCE_SQ;
    static double r_ext;
    
    static double P_MIN;
    static double P_MAX;
    
};

#endif /* PACKER_H */

