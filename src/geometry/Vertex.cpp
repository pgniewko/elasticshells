#include "Vertex.h"

//Vertex::Vertex() 
//{
//}

Vertex::Vertex(double x, double y, double z) : xyz(x, y, z)
{
//    velocity(0.0,0.0,0.0);
//    force(0.0,0.0,0.0);
    
    bondedNo = 0;
    noTris = 0;
    mass = 1.0;
    id = -1;
    for (int i = 0; i < NEIGH_SIZE; i++)
    {
        neigsIds[i] = -1.0;
        k0[i] = -1.0;
        trisIds[i] = -1;
    }
}

Vertex::Vertex(Vector3D v) : xyz(v) {}

Vertex::Vertex(const Vertex& orig) {}

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

void Vertex::addNeigh(int i)
{
    neigsIds[bondedNo] = i;
    bondedNo++;
}

void Vertex::addTriangle(int i)
{
    trisIds[noTris] = i;
    noTris++;
}
