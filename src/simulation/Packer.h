#ifndef PACKER_H
#define PACKER_H

#include "Environment.h"
#include "Vector3D.h"
#include "Cell.h"
#include "simulation/Box.h"

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
};

class Packer {
public:
    Packer();
    Packer(const Packer& orig);
    virtual ~Packer();
    
    static void packCells(Box&, std::vector<Cell>&);
    
private:
    
    static void fire(point_t*, int);
    static void velocityVerlet(point_t*, int);
    static void calcForces(point_t*, box_t&, int);
    static bool check_min_force(point_t* , int);
    
    static int FIRE_Nmin;
    static int FIRE_N;
    static double FIRE_DT;
    static double FIRE_ALPHA;
    static double FIRE_DTMAX;
    
    static double MIN_FORCE_SQ;
    static double delta_r;
    
};

#endif /* PACKER_H */

