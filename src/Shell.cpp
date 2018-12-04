#include "Shell.h"

utils::Logger Shell::cell_log("shell");

bool Shell::no_bending = false;

double Shell::FORCE_FRAC(0.0);
double Shell::MIN_FORCE(0.0);

Shell::Shell() {}

Shell::Shell(std::list<Triangle> tris) : shell_id(-1), number_v(0), number_t(0), number_s(0), nRT(0),
    V0(0)
{
    Tinker::constructVertices(*this, tris);
    Tinker::constructVTriangles(*this, tris);
    Tinker::constructTopology(*this);
    Tinker::constructBSprings(*this);
    randomRotate();
}

Shell::Shell(const Shell& orig) : center_of_mass(orig.center_of_mass), vertices(orig.vertices), triangles(orig.triangles), bhinges(orig.bhinges),
    shell_id(orig.shell_id), params(orig.params), number_v(orig.number_v), number_t(orig.number_t), number_s(orig.number_s),
    nRT(orig.nRT), V0(orig.V0) {}

Shell::~Shell() {}

double Shell::calcSurfaceArea(double d_param) const
{
    double totalSurface = 0.0;

    for (int i = 0; i < number_t; i++)
    {
        totalSurface += triangles[i].area(vertices, center_of_mass, d_param);
    }

    return totalSurface;
}

double Shell::calcSurfaceArea() const
{
    double surface = 0.0;

    for (int i = 0; i < number_t; i++)
    {
        surface += triangles[i].area(vertices);
    }

    return surface;
}

double Shell::calcVolume(double eps) const
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
    Vector3D tmp_m(0.0, 0.0, 0.0);
    double Mm = 0.0;

    for (int i = 0; i < number_v; i++)
    {
        tmp_m += vertices[i].r_c;
        Mm += 1.0;
    }

    if (Mm > 0.0)
    {
        tmp_m /= Mm;
        center_of_mass = tmp_m;
    }
    else
    {
        // REPORT PROBLEM
    }
}

void Shell::setBSprings(double E, double t, double nu_)
{
    if ( number_v == 1 || number_t == 1)
    {
        return;
    }

    for (int i = 0; i < number_s; i++)
    {
        bhinges[i].setD(E, t, nu_);
        bhinges[i].setThetaZero(vertices);
    }
}

void Shell::add_vector(const Vector3D& nxyz)
{
    for (int i = 0; i < number_v; i++)
    {
        vertices[i].r_c += nxyz;
    }
}

int Shell::getNumberTriangles() const
{
    return number_t;
}

int Shell::getNumberVertices() const
{
    return number_v;
}

int Shell::getNumberHinges() const
{
    return number_s;
}

void Shell::setVertexR(double rv)
{
    params.vertex_r = rv;
}

void Shell::setEcc(double a)
{
    params.ecc = a;
}

void Shell::setNu(double nu)
{
    params.nu = nu;
}

void Shell::setDp(double dP)
{
    setDp(dP, 0.0);
}

void Shell::setDp(double dP, double ddp)
{
    double randu = uniform(-ddp, ddp);
    params.dp = dP + randu;
    V0 = calcVolume();
    nRT = params.dp * V0 * ( 1.0 - OsmoticForce::getEpsilon() );
}

void Shell::setSpringConst(double E, double t, double nu_, std::string model_t)
{
    for (int i = 0; i < number_t; i++)
    {
        triangles[i].setParams(vertices, E, nu_, t);
    }
}

void Shell::setShellId(int ix)
{
    shell_id = ix;

    for (int i = 0; i < number_v; i++)
    {
        vertices[i].set_shell_id(shell_id);
    }
}

void Shell::setInitR(double rinit)
{
    params.init_r = rinit;
}

double Shell::getInitR() const
{
    return params.init_r;
}

Vector3D Shell::getCm() const
{
    return center_of_mass;
}

double Shell::getVertexR() const
{
    return params.vertex_r;
}

double Shell::getE() const
{
    return params.ecc;
}

double Shell::getNu() const
{
    return params.nu;
}

void Shell::randomRotate()
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

double Shell::getTurgor() const
{
    double turgor;

    if ( OsmoticForce::getFlag() )
    {
        double excluded_volume = V0 * OsmoticForce::getEpsilon();
        double cellVolume = calcVolume();

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

void Shell::update(double d)
{
    calc_cm();
}

void Shell::setConstantVolume(double scale)
{
    params.vol_c = V0 * (scale * scale * scale);

    if (calcVolume() != V0)
    {
        cell_log << utils::LogLevel::WARNING << "(calcVolume() != V0) @setConstantVolume\n";
    }
}

double Shell::checkVolumeCondition()
{
    double V = calcVolume();
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

std::ostream& operator<< (std::ostream& out, const Shell& c)
{
    out << "SHELL " << c.shell_id << ' ';
    out << c.number_v << ' ' << c.number_t << ' ' << c.number_s << ' ';
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

    for (int i = 0; i < c.number_s; i++)
    {
        out << "SHELLHINGE " <<  c.shell_id << ' ' << c.bhinges[i].getId() << ' ' << c.bhinges[i] << '\n';
    }

    return out;
}