#include "BendingSpring.h"
#include "Vertex.h"

BendingSpring::BendingSpring() : x1(-1), x2(-1), x3(-1), x4(-1) {}

BendingSpring::BendingSpring(int x1_, int x2_, int x3_, int x4_) : x1(x1_), x2(x2_), x3(x3_), x4(x4_)
{
}

BendingSpring::BendingSpring(const BendingSpring& orig) : D(orig.D), sinHalfTheta0(orig.sinHalfTheta0), x1(orig.x1), x2(orig.x2), x3(orig.x3), x4(orig.x4)
{}

BendingSpring::~BendingSpring() 
{}

void BendingSpring::setD(const double& E, const double& t, const double& nu)
{
    D = E*t*t*t / (12.0 * (1.0 - nu*nu));
}

void BendingSpring::calcBendingForces(Vertex vs[]) const
{
    
    Vector3D E = vs[x4].r_c - vs[x3].r_c;
    double E_norm = E.length();
    double E_norm2 = E_norm * E_norm;
    
    Vector3D N1 = 0.5 * cross(vs[x1].r_c - vs[x3].r_c, vs[x1].r_c - vs[x4].r_c);
    double N1_norm = N1.length();
    Vector3D N2 = 0.5 * cross(vs[x2].r_c - vs[x4].r_c, vs[x2].r_c - vs[x3].r_c);
    double N2_norm = N2.length();
    
    double sinHalfTheta = calcSinTheta(vs);
    
    double C = D * E_norm2 / (N1_norm + N2_norm) * (sinHalfTheta - sinHalfTheta0);
    
    Vector3D u1 = E_norm * N1 / (N1_norm*N1_norm);
    Vector3D u2 = E_norm * N2 / (N2_norm*N2_norm);
    
    Vector3D u3 =  dot(vs[x1].r_c-vs[x4].r_c, E) * N1 / (N1_norm*N1_norm*E_norm) + dot(vs[x2].r_c-vs[x4].r_c, E) * N2 / (N2_norm*N2_norm*E_norm);
    Vector3D u4 = -dot(vs[x1].r_c-vs[x3].r_c, E) * N1 / (N1_norm*N1_norm*E_norm) - dot(vs[x2].r_c-vs[x3].r_c, E) * N2 / (N2_norm*N2_norm*E_norm);
    
    vs[x1].f_c += (C * u1);
    vs[x2].f_c += (C * u2);
    vs[x3].f_c += (C * u3);
    vs[x4].f_c += (C * u4);   
    
}
    
void BendingSpring::setThetaZero(const Vertex vs[])
{
    sinHalfTheta0 = calcSinTheta(vs);
    if (sinHalfTheta0 < 0) // ENFORCE THAT THE INDEXING IS SUCH THAT THE ANGLE IS PI-THETA; THETA > 0 
    {
        int x_tmp = x1;
        x1 = x2;
        x2 = x_tmp;
    }
    sinHalfTheta0 = calcSinTheta(vs);

}

double BendingSpring::calcSinTheta(const Vertex vs[]) const
{
    
    Vector3D E = vs[x4].r_c - vs[x3].r_c;
    Vector3D e = E / E.length();

    Vector3D N1 = cross(vs[x3].r_c - vs[x1].r_c, vs[x4].r_c - vs[x1].r_c);
    Vector3D n1 = N1 / N1.length();
    Vector3D N2 = cross(vs[x4].r_c - vs[x2].r_c, vs[x3].r_c - vs[x2].r_c);
    Vector3D n2 = N2 / N2.length();

    int sign = SIGN( dot(cross(n1,n2),e) );
    
    double in_v = std::max(0.0, 0.5 * (1.0 - dot(n1,n2)) );
    return ( sign * sqrt( in_v ) );
}
