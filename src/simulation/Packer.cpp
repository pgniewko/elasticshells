#include "Packer.h"

double Packer::MIN_FORCE_SQ(1.0e-12);
double Packer::delta_r(1.0e-2);

int Packer::FIRE_Nmin(5);
int Packer::FIRE_N(0);
double Packer::FIRE_DT(0.0);
double Packer::FIRE_ALPHA(0.1);
double Packer::FIRE_DTMAX(0.0);
    
Packer::Packer() {
}

Packer::Packer(const Packer& orig) {
}

Packer::~Packer() {
}

void Packer::packCells(Box& boc, std::vector<Cell>& cells)
{
    int n = cells.size();
    point_t* points = new point_t[n];
    
    do
    {
        
        do
        {
            fire(points, n);
        }
        while ( check_min_force() );
    }
    while(true); // warunek jammingu, niezerowe cisnienie, bardzo male
    
    delete points;
    
}

void Packer::fire(point_t* points, box_t& box, int n)
{
    double f_inc = 1.1;
    double f_dec = 0.5;
    double a_start = 0.1;
    double f_a = 0.99;
    
    // MD step
    velocityVerlet(points, box, n); 
   
        // CALC P PARAMETER
    double P = 0.0;
    for (int i = 0; i < n; i++)
    {
        P += dot( points[i].f_c, points[i].v_c);
    }
    
    //=========================
    
    for (int i = 0; i < n; i++)
    {

        double v_norm = points[i].v_c.length();
        Vector3D F = points[i].f_c;
        F.normalize();
            
        points[i].v_c *= (1 - FIRE_ALPHA);
        points[i].v_c += FIRE_ALPHA * F * v_norm;
    }
    
    
    if (P > 0 && FIRE_N > FIRE_Nmin)
    {
        FIRE_DT = std::min(FIRE_DTMAX, FIRE_DT*f_inc);
        FIRE_ALPHA *= f_a;
    }
    
    FIRE_N++;
    
    if (P <= 0.0)
    {
        FIRE_DT *= f_dec;
        FIRE_ALPHA = a_start;
        
        for (int i = 0; i < n; i++)
        {
            points[i].v_c *= 0.0; // freeze the system
        }
        FIRE_N = 0;
    }
}

void Packer::velocityVerlet(point_t* points, box_t& box, int n)
{
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
    
    calcForces(points, box, n);
    // UPDATE VELOCITIES
    for (int i = 0; i < n; i++)
    {
        points[i].v_c += 0.5 * dt * points[i].f_c;
        points[i].f_p = points[i].f_c; // copy forces for the next time step integration
    }  
}

bool Packer::check_min_force(point_t* points, int n)
{
    for (int i = 0; i < n; i++)
    {

        if (points[i].f_c.length_sq() > MIN_FORCE_SQ)
        {
            return true;
        }

    }
    
    return false;
}