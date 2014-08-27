#include "VertexTriangle.h"

VertexTriangle::VertexTriangle()
{
}

VertexTriangle::VertexTriangle(int a, int b, int c) : ia(a), ib(b), ic(c) {}

VertexTriangle::VertexTriangle(const VertexTriangle& orig) : ia(orig.ia), ib(orig.ib), ic(orig.ic) 
{}

VertexTriangle::~VertexTriangle() {}

void VertexTriangle::setId(int idx)
{
    id = idx;
}

void VertexTriangle::printVertexTriangle()
{
    cout << "my id=" << id << " ";
    cout << " ia =" << ia << " ib =" << ib << " ic =" << ic << endl;
}

double VertexTriangle::area(const Vertex vs[])
{
    Triangle t(vs[ia].xyz, vs[ib].xyz, vs[ic].xyz);
    return t.area();
}