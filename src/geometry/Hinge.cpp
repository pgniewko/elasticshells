#include "Hinge.h"

Hinge::Hinge() : x1(-1), x2(-1), x3(-1), x4(-1), myid(-1) {}

Hinge::Hinge(int x1_, int x2_, int x3_, int x4_) : x1(x1_), x2(x2_), x3(x3_), x4(x4_), myid(-1)
{
}

Hinge::Hinge(const Hinge& orig) : D(orig.D),
    sinTheta0(orig.sinTheta0),
    theta0(orig.theta0),
    x1(orig.x1),
    x2(orig.x2),
    x3(orig.x3),
    x4(orig.x4),
    myid(orig.myid)
{
}

Hinge::~Hinge()
{
}

void Hinge::setD(const double& E, const double& t, const double& nu)
{
    D = E * t * t * t / (12.0 * (1.0 - nu * nu));
}

double Hinge::calcRadiusOfCurvature(std::vector<Vertex>& vs) const
{
    Vector3D E = vs[x4].r_c - vs[x3].r_c;
    double E_norm = E.length();

    Vector3D A1 = 0.5 * cross(vs[x1].r_c - vs[x3].r_c, vs[x1].r_c - vs[x4].r_c);
    double area1 = A1.length();

    Vector3D A2 = 0.5 * cross(vs[x2].r_c - vs[x4].r_c, vs[x2].r_c - vs[x3].r_c);
    double area2 = A2.length();

    double h1 = 2.0 * area1 / E_norm;
    double h2 = 2.0 * area2 / E_norm;
    double h = (h1 + h2) / 2.0;

    double R1 = 2.0 * fastmath::fast_sin(theta0 / 2.0) / h;
    R1 = 1.0 / R1;

    return R1;
}

void Hinge::setThetaZero(const std::vector<Vertex>& vs)
{
    sinTheta0 = calcSinTheta(vs);

    if (sinTheta0 < 0) // ENFORCE THAT THE INDEXING IS SUCH THAT THE ANGLE IS PI-THETA; THETA > 0
    {
        int x_tmp = x1;
        x1 = x2;
        x2 = x_tmp;
    }

    sinTheta0 = calcSinTheta(vs);
    theta0 = asin(sinTheta0);

}

double Hinge::calcTheta(const std::vector<Vertex>& vs) const
{
    double sin_theta = calcSinTheta(vs);
    return asin(sin_theta);
}

double Hinge::calcSinTheta(const std::vector<Vertex>& vs) const
{

    Vector3D E = vs[x4].r_c - vs[x3].r_c;
    Vector3D e = E / E.length();

    Vector3D N1 = cross(vs[x3].r_c - vs[x1].r_c, vs[x4].r_c - vs[x1].r_c);
    Vector3D n1 = N1 / N1.length();
    Vector3D N2 = cross(vs[x4].r_c - vs[x2].r_c, vs[x3].r_c - vs[x2].r_c);
    Vector3D n2 = N2 / N2.length();

    Vector3D cross_n1n2 = cross(n1, n2);
    int sign = SIGN( dot( cross_n1n2, e) );

    return ( sign * cross_n1n2.length() );
}

void Hinge::setId(int idx)
{
    myid = idx;
}

int Hinge::getId() const
{
    return myid;
}

std::ostream& operator<< (std::ostream& out, const Hinge& bs)
{
    out << bs.D << ' ' << bs.sinTheta0 << ' ' << bs.theta0 << ' ';
    out << bs.x1 << ' ' << bs.x2 << ' ' << bs.x3 << ' ' << bs.x4 << ' ';

    return out;
}