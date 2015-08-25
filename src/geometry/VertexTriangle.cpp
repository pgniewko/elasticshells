#include "VertexTriangle.h"

VertexTriangle::VertexTriangle() : ia(-1), ib(-1), ic(-1), myindex(-1) {}

VertexTriangle::VertexTriangle(int a, int b, int c) : ia(a), ib(b), ic(c), myindex(-1) {}

VertexTriangle::VertexTriangle(const VertexTriangle& orig) : ia(orig.ia), ib(orig.ib), ic(orig.ic), myindex(orig.myindex) {}

VertexTriangle::~VertexTriangle() {}

void VertexTriangle::setId(int idx)
{
    myindex = idx;
}

int VertexTriangle::getId()
{
    return myindex;
}

void VertexTriangle::printVertexTriangle()
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

double VertexTriangle::area(const Vertex vs[])
{
    if (ia != -1 && ib != -1 && ic != -1)
    {
        Triangle t(vs[ia].xyz, vs[ib].xyz, vs[ic].xyz);
        return t.area();
    }
    else
    {
        return 0.0;
    }
}

double VertexTriangle::area(const Vertex vs[], const Vector3D cm, double eps)
{
    if (ia != -1 && ib != -1 && ic != -1)
    {
        Vector3D ta = vs[ia].xyz - cm;
        Vector3D tb = vs[ib].xyz - cm;
        Vector3D tc = vs[ic].xyz - cm;
        double n_ta = ta.length() + eps;
        double n_tb = tb.length() + eps;
        double n_tc = tc.length() + eps;
        ta.setLength(n_ta);
        tb.setLength(n_tb);
        tc.setLength(n_tc);
        Triangle t(ta, tb, tc);
        return t.area();
    }
    else
    {
        return 0.0;
    }
}

Vector3D VertexTriangle::normal(const Vertex vs[])
{
    Vector3D ta = vs[ia].xyz;
    Vector3D tb = vs[ib].xyz;
    Vector3D tc = vs[ic].xyz;
    Triangle t(ta, tb, tc);
    Vector3D normal = t.normal();
    return normal;
}