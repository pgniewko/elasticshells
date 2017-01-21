#include "Packer.h"

double Packer::MIN_FORCE_SQ(1.0e-10);
double Packer::r_ext(1.0e-2);
double Packer::P_MIN(1e-5);
double Packer::P_MAX(2e-5);

int Packer::FIRE_Nmin(5);
int Packer::FIRE_N(0);
double Packer::FIRE_DT(0.1);
double Packer::FIRE_ALPHA(0.1);
double Packer::FIRE_DTMAX(1.0);
double Packer::FIRE_ITERS(0);
    
Packer::Packer() {
}

Packer::Packer(const Packer& orig) {}

Packer::~Packer() {
}

void Packer::packCells(Box& box, std::vector<Cell>& cells, double thickness)
{
    int n = cells.size();
    std::vector<point_t> points;
    for (int i = 0; i < n; i++)
    {
        point_t new_point;
        points.push_back(new_point);
    }
    
    box_t sim_box;
    sim_box.x = box.getXmin();
    sim_box.y = box.getYmin();
    sim_box.z = box.getZmin();
    sim_box.pbc = box.pbc;
    
    double radius_i;
    
    double E, t, P, r0, rv, nu;
    
    for (int i = 0; i < n; i++)
    {
        points[i].r_c.x = uniform(-sim_box.x, sim_box.x);
        points[i].r_c.y = uniform(-sim_box.y, sim_box.y);
        points[i].r_c.z = uniform(-sim_box.z, sim_box.z);
    }
    
    for (int i = 0; i < n; i++)
    {
        E = cells[i].getE();
        t = thickness;
        P = cells[i].getTurgor();
        r0 = cells[i].getInitR();
        rv = cells[i].getVertexR();
        nu = cells[i].getNu();
        radius_i = r0 * ( 1 + (1-nu) * P * r0 / (2.0*E*t)) + rv;
        points[i].radius = 0.01 * radius_i;
        points[i].radius_f = radius_i;
    }
    
    int loop_counter = 1;
    //return;
    
    
    do
    {
        std::cout << "-------------\nLOOP NUMBER=" << loop_counter << std::endl;
        Packer::inflatePoints(points);
        Packer::recenterCells(points, sim_box);

        
        std::cout << "RUNNING FIRE" << std::endl;
        do
        {
            Packer::fire(points, sim_box);
        }
        while ( Packer::check_min_force(points, sim_box) );
       
        
        //std::cout <<points.size()  << std::endl;
        
        for(int i = 0; i < points.size(); i++)
        {
            std::cout <<   points[i].r_c.x << " " << points[i].r_c.y << " " << points[i].r_c.z << " " << points[i].radius<< std::endl;
            //std::cout <<  "H "<<points[i].r_c.x << " " << points[i].r_c.y << " " << points[i].r_c.z << std::endl;
        }
        
        Packer::FIRE_ITERS = 0;
        Packer::FIRE_DT = 0.1;
        Packer::FIRE_ALPHA = 0.1;
        Packer::FIRE_N = 0;
        
        loop_counter++;
        
        for (int i = 0; i < n; i++)
        {
            points[i].v_c *= 0.0; // freeze the system
        }

        if ( points[0].radius > 2.0 || loop_counter > 400)
            break;
       
        //std::cout << "----------------------" << loop_counter << std::endl;
    }
    while( !Packer::jammed(points, sim_box) ); // warunek jammingu, niezerowe cisnienie, bardzo male

    
//    double box_scale = points[0].radius_f / points[0].radius;
//    
//    box.setX(box_scale * sim_box.x);
//    box.setY(box_scale * sim_box.y);
//    box.setZ(box_scale * sim_box.z);
//    for (int i = 0; i < points.size(); i++)
//    {
//        //assert( box_scale == points[i].radius_f / points[i].radius);
//        std::cout << box_scale<< " i: " << i << "r_f/r=" << points[i].radius_f / points[i].radius << std::endl;
//    }
}

void Packer::fire(std::vector<point_t>& points, box_t& box)
{
    Packer::FIRE_ITERS++;
    std::cout << "Packer::FIRE_DT=" << Packer::FIRE_DT << " Packer::FIRE_N=" << Packer::FIRE_N << " Packer::FIRE_ALPHA=" << Packer::FIRE_ALPHA  << std::endl;
//    //std::cout << "FIRE "<< std::endl;
    int n = points.size();
    double f_inc = 1.1;
    double f_dec = 0.5;
    double a_start = 0.1;
    double f_a = 0.99;
    
    // MD step
    Packer::velocityVerlet(points, box); 
//   
//        // CALC P PARAMETER
    double P = 0.0;
    for (int i = 0; i < n; i++)
    {
        P += dot( points[i].f_c, points[i].v_c);
        std::cout << " i="<< i << " P=" << P << " f_c=" << points[i].f_c << " v_c=" <<  points[i].v_c;
    }
    std::cout << std::endl;
    
    //=========================
    
    for (int i = 0; i < n; i++)
    {

        double v_length = points[i].v_c.length();
        Vector3D F = points[i].f_c;
        //std::cout << "BEFORE F=" << F << std::endl;
        
        F.normalize();
        
        //std::cout << "AFTER F=" << F << std::endl;
            
        points[i].v_c *= (1 - Packer::FIRE_ALPHA);
        points[i].v_c += Packer::FIRE_ALPHA * F * v_length;
    }
    
    
    if (P > 0 && Packer::FIRE_N > Packer::FIRE_Nmin)
    {
        Packer::FIRE_DT = std::min(Packer::FIRE_DTMAX, Packer::FIRE_DT*f_inc);
        Packer::FIRE_ALPHA *= f_a;
    }
    
    Packer::FIRE_N++;
    
    if (P <= 0.0)
    {
        
        Packer::FIRE_DT *= f_dec;
        Packer::FIRE_ALPHA = a_start;
        //std::cout <<"FIRE_DT=" << Packer::FIRE_DT  << std::endl;
        for (int i = 0; i < n; i++)
        {
            points[i].v_c *= 0.0; // freeze the system
        }
        Packer::FIRE_N = 0;
    }
}

void Packer::getDistance(Vector3D& dij, const Vector3D& vi, const Vector3D& vj, const box_t& box)
{
    dij = vj - vi;

    if (box.pbc)
    {
        double x, y, z;
        double bsx = 2 * box.x;
        double bsy = 2 * box.y;
        double bsz = 2 * box.z;
        x = round(dij.x / bsx) *  bsx;
        y = round(dij.y / bsy) *  bsy;
        z = round(dij.z / bsz) *  bsz;
        dij.x -= x;
        dij.y -= y;
        dij.z -= z;
    }
}

void Packer::calcBoxForces(std::vector<point_t>& points, const box_t& box)
{
    int n = points.size();
    Vector3D wallYZ, wallXZ, wallXY;
    Vector3D dij;
    double sgnx, sgny, sgnz;
    double bsx = box.x;
    double bsy = box.y;
    double bsz = box.z;
    double eb  = 1.0;
    double nub = 0.5;
    double rb_ = 0.0;
    double e1 = 1.0;
    double nu1 = 0.5;

    double r1;
    
    for (int i = 0; i < n; i++)
    {
        r1 = points[i].radius;
        sgnx = SIGN(points[i].r_c.x);
        wallYZ.x = sgnx * bsx;
        wallYZ.y = points[i].r_c.y;
        wallYZ.z = points[i].r_c.z;
        dij = points[i].r_c - wallYZ;
        points[i].f_c += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);

        sgny = SIGN(points[i].r_c.y);
        wallXZ.x = points[i].r_c.x;
        wallXZ.y = sgny * bsy;
        wallXZ.z = points[i].r_c.z;
        dij = points[i].r_c - wallXZ;
        points[i].f_c += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);

        sgnz = SIGN(points[i].r_c.z);
        wallXY.x = points[i].r_c.x;
        wallXY.y = points[i].r_c.y;
        wallXY.z = sgnz * bsz;
        dij = points[i].r_c - wallXY;
        points[i].f_c += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
        
        
        sgnx = SIGN(points[i].r_c.x);
        wallYZ.x = -sgnx * bsx;
        wallYZ.y = points[i].r_c.y;
        wallYZ.z = points[i].r_c.z;
        dij = points[i].r_c - wallYZ;
        points[i].f_c += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);

        sgny = SIGN(points[i].r_c.y);
        wallXZ.x = points[i].r_c.x;
        wallXZ.y = -sgny * bsy;
        wallXZ.z = points[i].r_c.z;
        dij = points[i].r_c - wallXZ;
        points[i].f_c += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);

        sgnz = SIGN(points[i].r_c.z);
        wallXY.x = points[i].r_c.x;
        wallXY.y = points[i].r_c.y;
        wallXY.z = -sgnz * bsz;
        dij = points[i].r_c - wallXY;
        points[i].f_c += HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
        
    }
}

void Packer::calcForces(std::vector<point_t>& points, box_t& box)
{
    int n = points.size();
    
    for (int i = 0; i < n; i++)
    {
        points[i].f_c.x = 0.0;
        points[i].f_c.y = 0.0;
        points[i].f_c.z = 0.0;
    }
    
    Vector3D force_ij;
    Vector3D dij;
    for (int i = 0; i < n; i++)
    {
        for (int j = i+1; j < n; j++)
        {
            
            Packer::getDistance(dij, points[j].r_c, points[i].r_c, box);
            force_ij = HertzianRepulsion::calcForce(dij, points[i].radius, points[j].radius, 1.0, 1.0, 0.5, 0.5);
            
            std::cout << "force_ij="<<force_ij << " distance="<< dij.length() << std::endl;
            points[i].f_c +=  force_ij;
            points[j].f_c += -force_ij; 
        }
    }
    
    if (!box.pbc)
    {
        calcBoxForces(points, box);
    }
}

void Packer::velocityVerlet(std::vector<point_t>& points, box_t& box)
{
    int n = points.size();
     double dt = FIRE_DT;

    // UPDATE POSITIONS
    for (int i = 0; i < n; i++)
    {
        points[i].r_c += dt * points[i].v_c + 0.5 * dt * dt * points[i].f_p; // use previously calculated forces
    }
    
    // UPDATE VELOCITIES
    for (int i = 0; i < n; i++)
    {

        points[i].v_c += 0.5 * dt * points[i].f_c;

    }
    
    calcForces(points, box);
    // UPDATE VELOCITIES
    for (int i = 0; i < n; i++)
    {

        points[i].v_c += 0.5 * dt * points[i].f_c;
        points[i].f_p = points[i].f_c; // copy forces for the next time step integration

    } 
}

bool Packer::check_min_force(std::vector<point_t>& points, box_t& box)
{
//    if (Packer::FIRE_ITERS > 1e3)
//    {
//        std::cout << "dupa" << std::endl;
//        exit(1);
//    }
        
    Packer::calcForces(points, box);
    
    int n = points.size();
    for (int i = 0; i < n; i++)
    {
        
        if (points[i].f_c.length_sq() > Packer::MIN_FORCE_SQ)
        {
            return true;
        }
    }
    
    return false;
}

bool Packer::jammed(std::vector<point_t>& points, box_t& box)
{
    int n = points.size();
    
    double pressure = Packer::calcPressure(points, box);
    
    std::cout << "jammed(...)" <<  Packer::P_MIN <<" " << pressure << " " << Packer::P_MAX << std::endl;
    
    if (pressure > Packer::P_MIN && pressure < Packer::P_MAX)
    {
        return true;
    }
    
    if (pressure > Packer::P_MAX && Packer::r_ext > 0.0 )
    {
        Packer::r_ext = -0.5 * Packer::r_ext;
    }
    
    if (pressure < Packer::P_MIN && Packer::r_ext < 0.0)
    {
        Packer::r_ext = -0.5 * Packer::r_ext;
    }
    
    std::cout << " New Packer::r_ext="<< Packer::r_ext <<std::endl;
    
    return false;
}

void Packer::inflatePoints(std::vector<point_t>& points)
{
    int n = points.size();
    
    for(int i = 0; i < n; i++)
    {
        std::cout << "inflatePoints =; i=" << i << " " <<points[i].radius << " ";
        points[i].radius *= (1.0 + Packer::r_ext);
        std::cout << points[i].radius << " Packer::r_ext="<< Packer::r_ext <<std::endl;
    }
}

void Packer::recenterCells(std::vector<point_t>& points, box_t& box)
{
    if (box.pbc)
    {
        double x, y, z;
        double bsx = 2 * box.x;
        double bsy = 2 * box.y;
        double bsz = 2 * box.z;
        for(uint i = 0; i < points.size(); i++)
        {
            x = round(points[i].r_c.x / bsx) *  bsx;
            y = round(points[i].r_c.y / bsy) *  bsy;
            z = round(points[i].r_c.z / bsz) *  bsz;
            Vector3D p_shift(-x, -y, -z);
            points[i].r_c += p_shift;
        }
    }
    
    return;
}

double Packer::calcPressure(std::vector<point_t>& points, box_t& box)
{
    double pressure;
    
    if (box.pbc)
    {
        double volume = 2.0 * box.x  * 2.0 * box.y * 2.0 * box.z;
        Vector3D rij;
        Vector3D fij;
        
        int n = points.size();
        for (int i = 0; i < n; i++)
        {
            for (int j = i+1; j < n; j++)
            {
                Packer::getDistance(rij, points[i].r_c, points[j].r_c, box);
                fij = HertzianRepulsion::calcForce(rij, points[i].radius, points[j].radius, 1.0, 1.0, 0.5, 0.5);
                pressure += dot(rij, fij);
            }
        }
        
        pressure /= (3*volume);
        return pressure;
    }
    else
    {
        double boxArea = 2.0 * box.x  * 2.0 * box.y * 2.0 * box.z;
        double totalForce = 0.0;
        
        for (uint i = 0; i < points.size(); i++)
        {
            totalForce += Packer::boxForce(points[i], box);
        }
        
        pressure = totalForce / boxArea; 
    }
    
    std::cout << "pressure=" << pressure<< std::endl;
    return pressure;
}

double Packer::boxForce(point_t& point, box_t& box)
{
    Vector3D wallYZ, wallXZ, wallXY;
    Vector3D dij;
    double sgnx, sgny, sgnz;
    double bsx = box.x;
    double bsy = box.y;
    double bsz = box.z;
    double eb  = 1.0;
    double nub = 0.5;
    double rb_ = 0.0;
    double e1 = 1.0;
    double nu1 = 0.5;
    
    Vector3D force;
    double force_collector = 0.0;
    
    double r1 = point.radius;
    sgnx = SIGN(point.r_c.x);
    wallYZ.x = sgnx * bsx;
    wallYZ.y = point.r_c.y;
    wallYZ.z = point.r_c.z;
    dij = point.r_c - wallYZ;
    force = HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
    force_collector += force.length();
    sgny = SIGN(point.r_c.y);
    wallXZ.x = point.r_c.x;
    wallXZ.y = sgny * bsy;
    wallXZ.z = point.r_c.z;
    dij = point.r_c - wallXZ;
    force = HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
    force_collector += force.length();
    sgnz = SIGN(point.r_c.z);
    wallXY.x = point.r_c.x;
    wallXY.y = point.r_c.y;
    wallXY.z = sgnz * bsz;
    dij = point.r_c - wallXY;
    force = HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
    force_collector += force.length();
    
    sgnx = SIGN(point.r_c.x);
    wallYZ.x = -sgnx * bsx;
    wallYZ.y = point.r_c.y;
    wallYZ.z = point.r_c.z;
    dij = point.r_c - wallYZ;
    force = HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
    force_collector += force.length();
    sgny = SIGN(point.r_c.y);
    wallXZ.x = point.r_c.x;
    wallXZ.y = -sgny * bsy;
    wallXZ.z = point.r_c.z;
    dij = point.r_c - wallXZ;
    force = HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
    force_collector += force.length();
    sgnz = SIGN(point.r_c.z);
    wallXY.x = point.r_c.x;
    wallXY.y = point.r_c.y;
    wallXY.z = -sgnz * bsz;
    dij = point.r_c - wallXY;
    force = HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);
    force_collector += force.length();
        
    //std::cout << "force_collector=" << force_collector << std::endl;
    return force_collector;
}