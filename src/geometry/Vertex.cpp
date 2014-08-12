#include "Vertex.h"

Vertex::Vertex() : xyz(0, 0, 0), mass(1.0)
{}

Vertex::Vertex(double x, double y, double z) : xyz(x, y, z), mass(1.0)
{}

Vertex::Vertex(Vector3D v) : xyz(v), mass(1.0) {}

Vertex::Vertex(const Vertex& orig) : xyz(orig.xyz) , id(orig.getId()), mass(orig.getMass())
{}

Vertex::~Vertex() {}

int Vertex::setId(int idx)
{
    id = idx;
    return id;
}

int Vertex::getId()
{
    return id;
}

double Vertex::setMass(double m)
{
    mass = m;
    return mass;
}

double Vertex::getMass()
{
    return mass;
}

void Vertex::printVertex()
{
    cout << " Memory address = " << this;
    cout << " Id = " << id;    
    cout << " Mass = " << mass;
    cout << " " << xyz;
    cout << endl;
}
