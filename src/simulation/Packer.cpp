#include "Packer.h"


utils::Logger Packer::packer_logs("packer");


double Packer::MIN_FORCE(1.0e-15);
double Packer::r_ext(1.0e-2);
double Packer::P_MIN(1e-9);
double Packer::P_MAX(2e-9);

int Packer::FIRE_Nmin(5);
int Packer::FIRE_N(0);
double Packer::FIRE_DT(0.1);
double Packer::FIRE_ALPHA(0.1);
double Packer::FIRE_DTMAX(1.5);

Packer::Packer() {}

Packer::Packer(const Packer& orig) {}

Packer::~Packer() {}

void Packer::packShells(Box& box, std::vector<Shell>& cells, double thickness, bool flag)
{
    std::size_t n = cells.size();
    std::vector<point_t> points;

    for (std::size_t i = 0; i < n; i++)
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

    double Z = 0.0;
    double P_final = 0.0;

    do
    {
        Packer::r_ext = 1.0e-2;

        for (std::size_t i = 0; i < n; i++)
        {
            points[i].r_c.x = uniform(-sim_box.x, sim_box.x);
            points[i].r_c.y = uniform(-sim_box.y, sim_box.y);
            points[i].r_c.z = uniform(-sim_box.z, sim_box.z);
            points[i].f_c *= 0.0;
            points[i].f_p *= 0.0;
            points[i].radius *= 0.0;
            points[i].radius_f;

        }

        for (std::size_t i = 0; i < n; i++)
        {
            E = cells[i].get_E();
            t = thickness;
            P = cells[i].get_turgor();
            r0 = cells[i].get_r0();
            rv = cells[i].get_vertex_size();
            nu = cells[i].get_nu();

            if (flag)
            {
                radius_i = r0 * ( 1 + (1 - nu) * P * r0 / (2.0 * E * t)) + rv;
            }
            else
            {
                radius_i = r0;
            }

            points[i].radius = 0.01 * radius_i;
            points[i].radius_f = radius_i;
        }


        do
        {
            Packer::inflatePoints(points);
            Packer::recenterShells(points, sim_box);

            do
            {
                Packer::fire(points, sim_box);
            }
            while ( Packer::check_min_force(points, sim_box) );

            Packer::FIRE_DT = 0.001;
            Packer::FIRE_ALPHA = 0.1;
            Packer::FIRE_N = 0;

            for (std::size_t i = 0; i < n; i++)
            {
                points[i].v_c *= 0.0; // freeze the system
            }
        }
        while ( !Packer::jammed(points, sim_box, P_final) ); // jamming condition, very small, residual pressure

    }
    while ( anyRattlerOrCrowder(points, sim_box, Z) );


    packer_logs << utils::LogLevel::INFO << "Jammed packing generated @:";
    packer_logs << utils::LogLevel::INFO << " P_MIN=" << Packer::P_MIN << " <= P=" << P_final << " <= P_MAX=" << Packer::P_MAX ;
    packer_logs << utils::LogLevel::INFO << " and <Z>=" << Z << "\n";


    double box_scale = points[0].radius_f / points[0].radius;

    box.set_x(box_scale * sim_box.x);
    box.set_y(box_scale * sim_box.y);
    box.set_z(box_scale * sim_box.z);

    for (std::size_t i = 0; i < n; i++)
    {
        Vector3D new_position = box_scale * points[i].r_c;
        cells[i].add_vector( -cells[i].get_cm() );
        cells[i].add_vector( new_position );
        cells[i].update();
    }
}

void Packer::fire(std::vector<point_t>& points, box_t& box)
{
    std::size_t n = points.size();
    double f_inc = 1.1;
    double f_dec = 0.25;
    double a_start = 0.1;
    double f_a = 0.99;

    // MD step
    Packer::velocityVerlet(points, box);


    // CALC P PARAMETER
    double P = 0.0;

    for (std::size_t i = 0; i < n; i++)
    {
        P += dot( points[i].f_c, points[i].v_c);
    }

    for (std::size_t i = 0; i < n; i++)
    {

        double v_length = points[i].v_c.length();
        Vector3D F = points[i].f_c;

        F.normalize();

        points[i].v_c *= (1 - Packer::FIRE_ALPHA);
        points[i].v_c += Packer::FIRE_ALPHA * F * v_length;
    }


    if (P > 0 && Packer::FIRE_N > Packer::FIRE_Nmin)
    {
        Packer::FIRE_DT = std::min(Packer::FIRE_DTMAX, Packer::FIRE_DT * f_inc);
        Packer::FIRE_ALPHA *= f_a;
    }

    Packer::FIRE_N++;

    if (P <= 0.0)
    {
        Packer::FIRE_DT *= f_dec;
        Packer::FIRE_ALPHA = a_start;

        for (std::size_t i = 0; i < n; i++)
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
    std::size_t n = points.size();
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

    for (std::size_t i = 0; i < n; i++)
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
    std::size_t n = points.size();

    for (std::size_t i = 0; i < n; i++)
    {
        points[i].f_c.x = 0.0;
        points[i].f_c.y = 0.0;
        points[i].f_c.z = 0.0;
    }

    Vector3D force_ij;
    Vector3D dij;

    for (std::size_t i = 0; i < n; i++)
    {
        for (std::size_t j = i + 1; j < n; j++)
        {

            Packer::getDistance(dij, points[j].r_c, points[i].r_c, box);
            force_ij = HertzianRepulsion::calcForce(dij, points[i].radius, points[j].radius, 1.0, 1.0, 0.5, 0.5);
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
    std::size_t n = points.size();
    double dt = FIRE_DT;

    // UPDATE POSITIONS
    for (std::size_t i = 0; i < n; i++)
    {
        points[i].r_c += dt * points[i].v_c + 0.5 * dt * dt * points[i].f_p; // use previously calculated forces
    }

    // UPDATE VELOCITIES
    for (std::size_t i = 0; i < n; i++)
    {
        points[i].v_c += 0.5 * dt * points[i].f_c;
    }

    calcForces(points, box);

    // UPDATE VELOCITIES
    for (std::size_t i = 0; i < n; i++)
    {

        points[i].v_c += 0.5 * dt * points[i].f_c;
        points[i].f_p = points[i].f_c; // copy forces for the next time step integration

    }
}

bool Packer::check_min_force(std::vector<point_t>& points, box_t& box)
{
    Packer::calcForces(points, box);
    std::size_t n = points.size();

    for (std::size_t i = 0; i < n; i++)
    {
        if (points[i].f_c.length() > Packer::MIN_FORCE)
        {
            return true;
        }
    }

    return false;
}

bool Packer::jammed(std::vector<point_t>& points, box_t& box, double& pf)
{
    double pressure = Packer::calcPressure(points, box);

    if (pressure > Packer::P_MIN && pressure < Packer::P_MAX)
    {
        pf = pressure;
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

    return false;
}

void Packer::inflatePoints(std::vector<point_t>& points)
{
    std::size_t n = points.size();

    for (std::size_t i = 0; i < n; i++)
    {
        points[i].radius *= (1.0 + Packer::r_ext);
    }
}

void Packer::recenterShells(std::vector<point_t>& points, box_t& box)
{
    if (box.pbc)
    {
        double x, y, z;
        double bsx = 2 * box.x;
        double bsy = 2 * box.y;
        double bsz = 2 * box.z;

        for (uint i = 0; i < points.size(); i++)
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
    double pressure = 0.0;

    if (box.pbc)
    {
        double volume = 2.0 * box.x  * 2.0 * box.y * 2.0 * box.z;
        Vector3D rij;
        Vector3D fij;

        std::size_t n = points.size();

        for (std::size_t i = 0; i < n; i++)
        {
            for (std::size_t j = i + 1; j < n; j++)
            {
                Packer::getDistance(rij, points[i].r_c, points[j].r_c, box);
                fij = HertzianRepulsion::calcForce(rij, points[i].radius, points[j].radius, 1.0, 1.0, 0.5, 0.5);
                pressure += dot(rij, fij);
            }
        }

        pressure /= (3.0 * volume);
        return pressure;
    }
    else
    {
        double boxArea = 2.0 * box.x  * 2.0 * box.y * 2.0 * box.z;
        double totalForce = 0.0;

        for (std::size_t i = 0; i < points.size(); i++)
        {
            totalForce += Packer::boxForce(points[i], box);
        }

        pressure = totalForce / boxArea;
    }

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

    return force_collector;
}

bool Packer::anyRattlerOrCrowder(const std::vector<point_t>& points, const box_t& box, double& Z)
{
    bool thereIsRattler = false;

    double contacts_sum = 0;


    std::size_t n = points.size();
    int* num_contacts = new int[n];

    for (std::size_t i = 0; i < n ; i++)
    {
        num_contacts[i] = 0;
    }

    for (std::size_t i = 0; i < n ; i++)
    {
        for (std::size_t j = 0; j < n; j++)
        {
            if (i != j)
            {
                num_contacts[i] += shellContacts(points[i], points[j], box);
            }
        }

        if (!box.pbc)
        {
            num_contacts[i] += boxContacts(points[i], box);
        }
    }

    for (std::size_t i = 0; i < n ; i++)
    {
        if (num_contacts[i] <= 3)  // simple criteria for rattlers
        {
            thereIsRattler = true;
        }

        contacts_sum += (double) num_contacts[i];
    }

    delete[] num_contacts;

    Z = contacts_sum / (double)n;

    double overpacked = false;

    if (Z >= 6.5 )
    {
        overpacked = true;
    }

    return (thereIsRattler || overpacked );
}

int Packer::boxContacts(const point_t& point, const box_t& box)
{
    Vector3D wallYZ, wallXZ, wallXY;
    Vector3D dij;
    int number_of_contacs = 0;
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

    double r1 = point.radius;
    sgnx = SIGN(point.r_c.x);
    wallYZ.x = sgnx * bsx;
    wallYZ.y = point.r_c.y;
    wallYZ.z = point.r_c.z;
    dij = point.r_c - wallYZ;
    force = HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);

    if ( force.length() )
    {
        number_of_contacs++;
    }


    sgny = SIGN(point.r_c.y);
    wallXZ.x = point.r_c.x;
    wallXZ.y = sgny * bsy;
    wallXZ.z = point.r_c.z;
    dij = point.r_c - wallXZ;
    force = HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);

    if ( force.length() )
    {
        number_of_contacs++;
    }

    sgnz = SIGN(point.r_c.z);
    wallXY.x = point.r_c.x;
    wallXY.y = point.r_c.y;
    wallXY.z = sgnz * bsz;
    dij = point.r_c - wallXY;
    force = HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);

    if ( force.length() )
    {
        number_of_contacs++;
    }

    sgnx = SIGN(point.r_c.x);
    wallYZ.x = -sgnx * bsx;
    wallYZ.y = point.r_c.y;
    wallYZ.z = point.r_c.z;
    dij = point.r_c - wallYZ;
    force = HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);

    if ( force.length() )
    {
        number_of_contacs++;
    }

    sgny = SIGN(point.r_c.y);
    wallXZ.x = point.r_c.x;
    wallXZ.y = -sgny * bsy;
    wallXZ.z = point.r_c.z;
    dij = point.r_c - wallXZ;
    force = HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);

    if ( force.length() )
    {
        number_of_contacs++;
    }

    sgnz = SIGN(point.r_c.z);
    wallXY.x = point.r_c.x;
    wallXY.y = point.r_c.y;
    wallXY.z = -sgnz * bsz;
    dij = point.r_c - wallXY;
    force = HertzianRepulsion::calcForce(dij, r1, rb_, e1, eb, nu1, nub);

    if ( force.length() )
    {
        number_of_contacs++;
    }


    return number_of_contacs;
}

int Packer::shellContacts(const point_t& point_1, const point_t& point_2, const box_t& box)
{
    Vector3D force_ij;
    Vector3D dij;

    Packer::getDistance(dij, point_2.r_c, point_1.r_c, box);
    force_ij = HertzianRepulsion::calcForce(dij, point_1.radius, point_2.radius, 1.0, 1.0, 0.5, 0.5);

    if ( force_ij.length_sq() > 0)
    {
        return 1;
    }

    return 0;
}