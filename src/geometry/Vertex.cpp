#include "Vertex.h"

Vertex::Vertex() : xyz(0, 0, 0), nneigh(0), ntrian(0), mass(1.0){}

Vertex::Vertex(double x, double y, double z) : xyz(x, y, z), nneigh(0), ntrian(0), mass(1.0){}

Vertex::Vertex(Vector3D v) : xyz(v), nneigh(0), ntrian(0), mass(1.0) {}

Vertex::Vertex(const Vertex& orig) : xyz(orig.xyz), nneigh(orig.nneigh), ntrian(orig.ntrian), id(orig.getId()), mass(orig.getMass())
{}

Vertex::~Vertex() {}

void Vertex::addNeighbor(int idx, double k0n)
{
    
    for (int i = 0; i < nneigh; i++)
    {
        if (neighbors[i] == idx)  return;
    }
    neighbors[nneigh] = idx;
    R0[nneigh] = k0n;
    nneigh++;
}

void Vertex::addTriangle(int idx)
{
    for (int i = 0; i < ntrian; i++)
    {
        if (vertextri[i] == idx) return;
    }
    vertextri[ntrian] = idx;
    ntrian++;    
}

bool Vertex::isBonded(int vidx)
{
    for (int i = 0; i < nneigh; i++)
    {
        if (neighbors[i] == vidx) return true;
    }
    return false;
}

void Vertex::voidForce()
{
    force.x = 0.0;
    force.y = 0.0;
    force.z = 0.0;
}

void Vertex::voidVelocity()
{
    velocity.x = 0.0;
    velocity.y = 0.0;
    velocity.z = 0.0;
}

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
    cout << "nneigh= " << nneigh << " ntrian=" << ntrian << " ";
    for (int i = 0; i < nneigh; i++)
    {
        cout << neighbors[i] << " " ;
    }
    cout << " : ";
    
    for (int i = 0; i < ntrian; i++)
    {
        cout << vertextri[i] << " " ;
    }   
    cout <<endl;
}
