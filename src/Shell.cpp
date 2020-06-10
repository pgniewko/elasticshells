#include "Shell.h"

utils::Logger Shell::cell_log("shell");

bool Shell::bending = false;

double Shell::FORCE_FRAC(0.0);
double Shell::MIN_FORCE(0.0);

Shell::Shell() {}

Shell::Shell(int nv, int nt, int nh)
{
    for (int i = 0; i < nv; i++)
    {
        Vertex new_vertx;
        vertices.push_back( new_vertx );
    }

    for (int i = 0; i < nt; i++)
    {
        Element new_vertx_triangle;
        triangles.push_back( new_vertx_triangle );
    }

    for (int i = 0; i < nh; i++)
    {
        Hinge new_hinge;
        hinges.push_back( new_hinge );
    }

    number_v = nv;
    number_t = nt;
    number_h = nh;
}

Shell::Shell(std::list<Triangle> tris) : shell_id(-1),
    number_v(0),
    number_t(0),
    number_h(0),
    nRT(0),
    V0(0)
{
    Tinker::construct_vertices(*this, tris);
    Tinker::construct_elements(*this, tris);
    Tinker::construct_topology(*this);
    Tinker::construct_hinges(*this);
    random_rotate();
}

Shell::Shell(const Shell& orig) : center_of_mass(orig.center_of_mass),
    vertices(orig.vertices),
    triangles(orig.triangles),
    hinges(orig.hinges),
    shell_id(orig.shell_id),
    params(orig.params),
    number_v(orig.number_v),
    number_t(orig.number_t),
    number_h(orig.number_h),
    nRT(orig.nRT), V0(orig.V0)
{
    for (int i = 0; i < number_v; i++)
    {
        vertices.push_back( orig.vertices[i] );
    }

    for (int i = 0; i < number_t; i++)
    {
        triangles.push_back( orig.triangles[i] );
    }

    for (int i = 0; i < number_h; i++)
    {
        hinges.push_back( orig.hinges[i] );
    }
}

Shell::~Shell()
{
}

double Shell::calc_volume(double eps) const
{
    double volume = 0.0;

    if (number_v == 1)
    {
        volume = 4.0 / 3.0 * constants::pi * params.init_r * params.init_r * params.init_r;
    }
    else
    {
        int va, vb, vc;

        for (int i = 0; i < number_t; i++)
        {
            va = triangles[i].ia;
            vb = triangles[i].ib;
            vc = triangles[i].ic;
            volume += Tetrahedron::volume(vertices[va].r_c, vertices[vb].r_c, vertices[vc].r_c, center_of_mass, eps);
        }
    }

    return volume;
}

void Shell::calc_cm()
{
    Vector3D tmp(0.0, 0.0, 0.0);
    double mass = 0.0;

    for (int i = 0; i < number_v; i++)
    {
        tmp += vertices[i].r_c;
        mass += 1.0;
    }

    if (mass > 0.0)
    {
        tmp /= mass;
        center_of_mass = tmp;
    }
    else
    {
        // REPORT PROBLEM
    }
}

void Shell::set_hinges(double E, double t, double nu_)
{
    if ( number_v == 1 || number_t == 1)
    {
        return;
    }

    for (int i = 0; i < number_h; i++)
    {
        hinges[i].set_d(E, t, nu_);
        hinges[i].set_theta_zero(vertices);
    }
}

void Shell::add_vector(const Vector3D& nxyz)
{
    for (int i = 0; i < number_v; i++)
    {
        vertices[i].r_c += nxyz;
    }
}

int Shell::get_number_triangles() const
{
    return number_t;
}

int Shell::get_number_vertices() const
{
    return number_v;
}

int Shell::get_number_hinges() const
{
    return number_h;
}

void Shell::set_vertex_size(double rv)
{
    params.vertex_r = rv;
}

void Shell::set_ecc(double a)
{
    params.ecc = a;
}

void Shell::set_nu(double nu)
{
    params.nu = nu;
}

void Shell::set_dp(double dP)
{
    set_dp(dP, 0.0);
}

void Shell::set_dp(double dP, double ddp)
{
    double randu = uniform(-ddp, ddp);
    params.dp = dP + randu;
    V0 = calc_volume();
    nRT = params.dp * V0 * (1.0 - OsmoticForce::get_epsilon());
}

void Shell::set_elements_parameters(double E, double t, double nu_)
{
    for (int i = 0; i < number_t; i++)
    {
        triangles[i].set_params(vertices, E, nu_, t);
    }
}

void Shell::set_shell_id(int ix)
{
    shell_id = ix;

    for (int i = 0; i < number_v; i++)
    {
        vertices[i].set_shell_id(shell_id);
    }
}

void Shell::set_r0(double rinit)
{
    params.init_r = rinit;
}

double Shell::get_r0() const
{
    return params.init_r;
}

Vector3D Shell::get_cm() const
{
    return center_of_mass;
}

double Shell::get_vertex_size() const
{
    return params.vertex_r;
}

double Shell::get_E() const
{
    return params.ecc;
}

double Shell::get_nu() const
{
    return params.nu;
}

double Shell::get_turgor() const
{
    double turgor;

    if ( OsmoticForce::get_flag() )
    {
        double excluded_volume = V0 * OsmoticForce::get_epsilon();
        double cellVolume = calc_volume();

        if ( (cellVolume - excluded_volume) > 0 )
        {
            turgor = nRT / (cellVolume - excluded_volume);
        }
        else
        {
            cell_log << utils::LogLevel::SEVERE << " CELL VOLEUME SMALLER THAN OSMOTIC EXCLUDED VOLUME !\n";
            turgor = nRT / (0.01 * V0) ;
        }
    }
    else
    {
        turgor = params.dp;
    }

    return turgor;
}

void Shell::set_constant_volume(double scale)
{
    params.vol_c = V0 * (scale * scale * scale);

    if (calc_volume() != V0)
    {
        cell_log << utils::LogLevel::WARNING << "(calcVolume() != V0) @setConstantVolume\n";
    }
}

double Shell::check_volume_condition()
{
    double V = calc_volume();
    return (params.vol_c - V) / V;
}

double Shell::ajust_turgor(double step)
{
    params.dp = (1.0 + step) * params.dp;
    return params.dp;
}

const shell_params_t& Shell::get_params() const
{
    return params;
}

double Shell::calc_surface_area(double d_param) const
{
    double totalSurface = 0.0;

    for (int i = 0; i < number_t; i++)
    {
        totalSurface += triangles[i].area(vertices, center_of_mass, d_param);
    }

    return totalSurface;
}

std::ostream& operator<< (std::ostream& out, const Shell& c)
{
    out << "SHELL " << c.shell_id << ' ';
    out << c.number_v << ' ' << c.number_t << ' ' << c.number_h << ' ';
    out << c.params.vertex_r << ' ' << c.params.ecc << ' ' << c.params.nu << ' ';
    out << c.params.dp << ' ' << c.params.init_r << ' ' << c.params.vol_c << ' ';
    out << c.nRT << ' ' << c.V0 << "\n";

    for (int i = 0; i < c.number_v; i++)
    {
        out << "SHELLVERTEX " <<  c.shell_id << ' ' << c.vertices[i] << '\n';
    }

    for (int i = 0; i < c.number_t; i++)
    {
        out << "SHELLTRIANG " <<  c.shell_id << ' ' << c.triangles[i] << '\n';
    }

    for (int i = 0; i < c.number_h; i++)
    {
        out << "SHELLHINGE " <<  c.shell_id << ' ' << c.hinges[i].get_id() << ' ' << c.hinges[i] << '\n';
    }

    return out;
}

void Shell::random_rotate()
{
    calc_cm();

    double u1 = uniform();
    double u2 = uniform();
    double u3 = uniform();
    double q0 = fastmath::fast_sqrt(1 - u1) * fastmath::fast_sin(2 * M_PI * u2);
    double q1 = fastmath::fast_sqrt(1 - u1) * fastmath::fast_cos(2 * M_PI * u2);
    double q2 = fastmath::fast_sqrt(u1) * fastmath::fast_sin(2 * M_PI * u3);
    double q3 = fastmath::fast_sqrt(u1) * fastmath::fast_cos(2 * M_PI * u3);
    double A[3][3];
    A[0][0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
    A[0][1] = 2 * (q1 * q2 + q0 * q3);
    A[0][2] = 2 * (q1 * q3 - q0 * q2);
    A[1][0] = 2 * (q1 * q2 - q0 * q3);
    A[1][1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
    A[1][2] = 2 * (q2 * q3 + q0 * q1);
    A[2][0] = 2 * (q1 * q3 + q0 * q2);
    A[2][1] = 2 * (q2 * q3 - q0 * q1);
    A[2][2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
    double xnew = 0.0;
    double ynew = 0.0;
    double znew = 0.0;

    for (int i = 0; i  < number_v; i++)
    {
        double xi = vertices[i].r_c.x - center_of_mass.x;
        double yi = vertices[i].r_c.y - center_of_mass.y;
        double zi = vertices[i].r_c.z - center_of_mass.z;
        xnew  = A[0][0] * xi;
        xnew += A[0][1] * yi;
        xnew += A[0][2] * zi;
        ynew  = A[1][0] * xi;
        ynew += A[1][1] * yi;
        ynew += A[1][2] * zi;
        znew  = A[2][0] * xi;
        znew += A[2][1] * yi;
        znew += A[2][2] * zi;
        vertices[i].r_c.x = xnew + center_of_mass.x;
        vertices[i].r_c.y = ynew + center_of_mass.y;
        vertices[i].r_c.z = znew + center_of_mass.z;
    }
}