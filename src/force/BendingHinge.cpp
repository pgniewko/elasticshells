#include "BendingHinge.h"

BendingHinge::BendingHinge() : x1(-1), x2(-1), x3(-1), x4(-1) 
{}

BendingHinge::BendingHinge(int x1_, int x2_, int x3_, int x4_) : x1(x1_), x2(x2_), x3(x3_), x4(x4_) 
{}

BendingHinge::BendingHinge(const BendingHinge& orig) : D(orig.D), sinTheta0(orig.sinTheta0), theta0(orig.theta0), x1(orig.x1), x2(orig.x2), x3(orig.x3), x4(orig.x4)
{}

BendingHinge::~BendingHinge()
{}

void BendingHinge::setD(const double& E, const double& t, const double& nu)
{
    D = E * t * t * t / (12.0 * (1.0 - nu * nu));
}

double BendingHinge::calcRadiusOfCurvature(Vertex vs[]) const
{
    Vector3D E = vs[x4].r_c - vs[x3].r_c;
    double E_norm = E.length();

    Vector3D A1 = 0.5 * cross(vs[x1].r_c - vs[x3].r_c, vs[x1].r_c - vs[x4].r_c);
    double area1 = A1.length();

    Vector3D A2 = 0.5 * cross(vs[x2].r_c - vs[x4].r_c, vs[x2].r_c - vs[x3].r_c);
    double area2 = A2.length();

    double h1 = 2.0 * area1 / E_norm;
    double h2 = 2.0 * area2 / E_norm;
    double h = (h1 + h2)/2.0;

    double R1 = 2.0 * fastmath::fast_sin(theta0 / 2.0) / h;
    R1 = 1.0 / R1;
    
    return R1;
}

void BendingHinge::calcBendingForces(Vertex vs[]) const
{
    
    Vector3D E = vs[x4].r_c - vs[x3].r_c;
    double E_norm = E.length();
    double E_norm2 = E_norm * E_norm;

    Vector3D A1 = 0.5 * cross(vs[x1].r_c - vs[x3].r_c, vs[x1].r_c - vs[x4].r_c);
    double area1 = A1.length();
    Vector3D n1 = A1 / area1;

    Vector3D A2 = 0.5 * cross(vs[x2].r_c - vs[x4].r_c, vs[x2].r_c - vs[x3].r_c);
    double area2 = A2.length();
    Vector3D n2 = A2 / area2;

    //double sinTheta = calcSinTheta(vs);
    double theta = calcTheta(vs);

    double C = D * E_norm2 / (area1 + area2) * fastmath::fast_sin(theta - theta0); // Grinspun eq.3, Cubic Shells, 2007


    Vector3D u1 = E_norm / (2.0 * area1) * n1;
    Vector3D u2 = E_norm / (2.0 * area2) * n2;

    Vector3D u3 =  dot(vs[x1].r_c - vs[x4].r_c, E) * n1 / (2.0 * area1 * E_norm) + dot(vs[x2].r_c - vs[x4].r_c, E) * n2 / (2.0 * area2 * E_norm);
    Vector3D u4 = -dot(vs[x1].r_c - vs[x3].r_c, E) * n1 / (2.0 * area1 * E_norm) - dot(vs[x2].r_c - vs[x3].r_c, E) * n2 / (2.0 * area2 * E_norm);

    vs[x1].f_c += (C * u1);
    vs[x2].f_c += (C * u2);
    vs[x3].f_c += (C * u3);
    vs[x4].f_c += (C * u4);

}

void BendingHinge::setThetaZero(const Vertex vs[])
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

double BendingHinge::calcTheta(const Vertex vs[]) const
{
    double sin_theta = calcSinTheta(vs);
    return asin(sin_theta);
}

double BendingHinge::calcSinTheta(const Vertex vs[]) const
{

    Vector3D E = vs[x4].r_c - vs[x3].r_c;
    Vector3D e = E / E.length();

    Vector3D N1 = cross(vs[x3].r_c - vs[x1].r_c, vs[x4].r_c - vs[x1].r_c);
    Vector3D n1 = N1 / N1.length();
    Vector3D N2 = cross(vs[x4].r_c - vs[x2].r_c, vs[x3].r_c - vs[x2].r_c);
    Vector3D n2 = N2 / N2.length();

    Vector3D cross_n1n2 = cross(n1, n2);
    int sign = SIGN( dot( cross_n1n2 , e) );

    return ( sign * cross_n1n2.length() );
}
