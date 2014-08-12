#include "VertexTriangle.h"

VertexTriangle::VertexTriangle()
{
}

VertexTriangle::VertexTriangle(Vertex * v, Vertex * w, Vertex * z) 
{
    a = v;
    b = w;
    c = z;
}

VertexTriangle::VertexTriangle(const VertexTriangle& orig) : a(orig.a), b(orig.b), c(orig.c) 
{
}

VertexTriangle::~VertexTriangle() {}

void VertexTriangle::setId(int idx)
{
    id = idx;
}

double VertexTriangle::area()
{
    Triangle t(a->xyz, b->xyz, c->xyz);
    return t.area();
}
