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
    cout << "my id=" << myindex << " ";
    cout << " ia =" << ia << " ib =" << ib << " ic =" << ic << endl;
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