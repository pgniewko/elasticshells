#include "VertexTriangle.h"

double cot(double x) {return 1.0 / tan(x);}

VertexTriangle::VertexTriangle() {}

VertexTriangle::VertexTriangle(int a, int b, int c) : ia(a), ib(b), ic(c), myindex(-1) {}

VertexTriangle::VertexTriangle(const VertexTriangle& orig) : ia(orig.ia), ib(orig.ib), ic(orig.ic), myindex(orig.myindex) 
{
    for (int i = 0; i < 3; i++)
    {
        an[i] = orig.an[i];
        L2[i] = orig.L2[i];
        ki[i] = orig.ki[i];
        ci[i] = orig.ci[i];
    }
}

VertexTriangle::~VertexTriangle() {}

void VertexTriangle::setId(int idx)
{
    myindex = idx;
}

int VertexTriangle::getId() const
{
    return myindex;
}

void VertexTriangle::printVertexTriangle() const
{
    std::cout << "my id=" << myindex << " ";
    std::cout << " ia =" << ia << " ib =" << ib << " ic =" << ic << std::endl;
}

void VertexTriangle::subsVertex(int ix_old, int ix_new)
{
    if (ix_old == ia)
    {
        ia = ix_new;
        return;
    }
    else if (ix_old == ib)
    {
        ib = ix_new;
        return;
    }
    else if (ix_old == ic)
    {
        ic = ix_new;
        return;
    }

    return;
}

double VertexTriangle::area(const Vertex vs[]) const
{
    if (ia != -1 && ib != -1 && ic != -1)
    {
        Triangle t(vs[ia].r_c, vs[ib].r_c, vs[ic].r_c);
        return t.area();
    }
    else
    {
        return 0.0;
    }
}

double VertexTriangle::area(const Vertex vs[], const Vector3D cm, double eps) const
{
    if (ia != -1 && ib != -1 && ic != -1)
    {
        Vector3D ta = vs[ia].r_c - cm;
        Vector3D tb = vs[ib].r_c - cm;
        Vector3D tc = vs[ic].r_c - cm;
        double n_ta = ta.length() + eps;
        double n_tb = tb.length() + eps;
        double n_tc = tc.length() + eps;
        ta.set_length(n_ta);
        tb.set_length(n_tb);
        tc.set_length(n_tc);
        Triangle t(ta, tb, tc);
        return t.area();
    }
    else
    {
        return 0.0;
    }
}

Vector3D VertexTriangle::normal(const Vertex vs[]) const
{
    Vector3D ta = vs[ia].r_c;
    Vector3D tb = vs[ib].r_c;
    Vector3D tc = vs[ic].r_c;
    Triangle t(ta, tb, tc);
    Vector3D normal = t.normal();
    return normal;
}

void VertexTriangle::setL2(const Vertex vs[])
{
    L2[0] = (vs[ib].r_c - vs[ic].r_c).length_sq();
    L2[1] = (vs[ia].r_c - vs[ic].r_c).length_sq();
    L2[2] = (vs[ia].r_c - vs[ib].r_c).length_sq();
}

void VertexTriangle::setAn(const Vertex vs[])
{
    // MAKE SURE THAT the angle is between 0-180
    Vector3D ca = vs[ic].r_c - vs[ia].r_c;
    Vector3D ba = vs[ib].r_c - vs[ia].r_c;
    an[0] = ca.angle(ba);
    
    Vector3D ab = vs[ia].r_c - vs[ib].r_c;
    Vector3D cb = vs[ic].r_c - vs[ib].r_c;
    an[1] = ab.angle(cb);

    Vector3D ac = vs[ia].r_c - vs[ic].r_c;
    Vector3D bc = vs[ib].r_c - vs[ic].r_c;
    an[2] = ac.angle(bc);
}

void VertexTriangle::setKi(const Vertex vs[], const double &E, const double &nu, const double &t)
{
    double Ap = area(vs);
    ki[0] = E * t * (2.0 * cot(an[0]) * cot(an[0]) + 1.0 - nu) / (16.0 * Ap * (1.0 - nu*nu));
    ki[1] = E * t * (2.0 * cot(an[1]) * cot(an[1]) + 1.0 - nu) / (16.0 * Ap * (1.0 - nu*nu));
    ki[2] = E * t * (2.0 * cot(an[2]) * cot(an[2]) + 1.0 - nu) / (16.0 * Ap * (1.0 - nu*nu));
    
}

void VertexTriangle::setCi(const Vertex vs[], const double &E, const double &nu, const double &t)
{
    double Ap = area(vs);
    ci[0] = E*t*(2.0*cot(an[1])*cot(an[2]) + nu - 1.0 ) / (16.0 * Ap * (1.0 - nu*nu));
    ci[1] = E*t*(2.0*cot(an[0])*cot(an[2]) + nu - 1.0 ) / (16.0 * Ap * (1.0 - nu*nu));
    ci[2] = E*t*(2.0*cot(an[0])*cot(an[1]) + nu - 1.0 ) / (16.0 * Ap * (1.0 - nu*nu));
}

void VertexTriangle::setParams(const Vertex vs[], const double E, const double nu, const double t)
{
    setL2(vs);
    setAn(vs);
    setKi(vs, E, nu, t);
    setCi(vs, E, nu, t);
}

void VertexTriangle::calcFemForces(Vertex vs[]) const
{
    // 1 - a; 2 - b; 3 - c;
    double l0_sq = (vs[ib].r_c - vs[ic].r_c).length_sq() - L2[0];
    double l1_sq = (vs[ia].r_c - vs[ic].r_c).length_sq() - L2[1];
    double l2_sq = (vs[ia].r_c - vs[ib].r_c).length_sq() - L2[2];
    
    Vector3D T11;
    Vector3D T12; 
    T11 += ki[2] * l2_sq * (vs[ib].r_c - vs[ia].r_c) + ki[1] * l1_sq * (vs[ic].r_c - vs[ia].r_c);
    T12 += (ci[1] * l0_sq + ci[0]*l1_sq) * (vs[ib].r_c - vs[ia].r_c);
    T12 += (ci[2] * l0_sq + ci[0]*l2_sq) * (vs[ic].r_c - vs[ia].r_c);
    
    vs[ia].f_c += (T11 + T12);
    
    Vector3D T21;
    Vector3D T22; 
    T21 += ki[2] * l2_sq * (vs[ia].r_c - vs[ib].r_c) + ki[0] * l0_sq * (vs[ic].r_c - vs[ib].r_c);
    T22 += (ci[0] * l1_sq + ci[1]*l0_sq) * (vs[ia].r_c - vs[ib].r_c);
    T22 += (ci[2] * l1_sq + ci[1]*l2_sq) * (vs[ic].r_c - vs[ib].r_c);
    
    vs[ib].f_c += (T21 + T22);
    
    Vector3D T31;
    Vector3D T32; 
    T31 += ki[1] * l1_sq * (vs[ia].r_c - vs[ic].r_c) + ki[0] * l0_sq * (vs[ib].r_c - vs[ic].r_c);
    T32 += (ci[0] * l2_sq + ci[2]*l0_sq) * (vs[ia].r_c - vs[ic].r_c);
    T32 += (ci[1] * l2_sq + ci[2]*l1_sq) * (vs[ib].r_c - vs[ic].r_c);
    
    vs[ic].f_c += (T31 + T32);
    
}

