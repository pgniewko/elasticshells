#include "Vertex.h"

Vertex::Vertex() : xyz(0, 0, 0), nneigh(0), ntrian(0), nbneigh(0), id(-1), mass(1.0) {}

Vertex::Vertex(double x, double y, double z) : xyz(x, y, z), nneigh(0), ntrian(0), nbneigh(0), id(-1), mass(1.0)
{}

Vertex::Vertex(Vector3D v) : xyz(v), nneigh(0), ntrian(0), nbneigh(0), id(-1), mass(1.0) {}

Vertex::Vertex(const Vertex& orig) : xyz(orig.xyz), force(orig.force), velocity(orig.velocity),
                                     tmp_xyz(orig.tmp_xyz), tmp_force(orig.tmp_force), tmp_velocity(orig.tmp_velocity),
                                     nneigh(orig.nneigh), ntrian(orig.ntrian), nbneigh(orig.nbneigh),
                                     id(orig.id), mass( orig.mass ), visc(orig.visc)
{
    for (int i = 0; i < nneigh; i++)
    {
        neighbors[i] = orig.neighbors[i];
        R0[i] = orig.R0[i];
    }
    
    for (int i = 0; i < ntrian; i++)
    {
        vertextri[i] = orig.vertextri[i];
    }
    
    for (int i = 0; i < nbneigh; i++)
    {
        nbvertices[i] = orig.nbvertices[i];
        nbcellid[i] = orig.nbcellid[i];
    }
}

Vertex::~Vertex() {}

void Vertex::addNeighbor(int idx, double k0n)
{
 
    try
    {
        if (getNumNeighbors() >= MAX_CELLS)
            throw MaxSizeException("Maximum number of neighbors has been reached."
                                    "New neighbor will not be added !");
        
        for (int i = 0; i < nneigh; i++)
        {
            if (neighbors[i] == idx)  return;
        }
        neighbors[nneigh] = idx;
        R0[nneigh] = k0n;
        nneigh++;
    }
    catch (MaxSizeException& e)
    {
        cout << e.what() << endl;
        return;
    }
}

bool Vertex::isNeighbor(int vidx)
{
    for (int i = 0; i < nneigh; i++)
    {
        if (neighbors[i] == vidx) return true;
    }
    return false;
}

void Vertex::addTriangle(int idx)
{
    try
    {
        if (getNumNeighbors() >= TRIAN_MAX)
            throw MaxSizeException("Maximum number of triangles has been reached."
                                    "New triangle will not be added !");
        
        for (int i = 0; i < ntrian; i++)
        {
            if (vertextri[i] == idx) return;
        }
        vertextri[ntrian] = idx;
        ntrian++;    
    }
    catch (MaxSizeException& e)
    {
        cout << e.what() << endl;
        return;
    }
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

int Vertex::getNumNeighbors()
{
    return nneigh;
}

int Vertex::getNumVTriangles()
{
    return ntrian;
}

int Vertex::getNeighborId(int idx)
{
    return neighbors[idx];
}

int Vertex::getTriangleId(int idx)
{
    return vertextri[idx];
}

double Vertex::getNeighborR0(int idx)
{
    return R0[idx];
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

double Vertex::setVisc(double v)
{
    visc = v;
    return visc;
}

double Vertex::getVisc()
{
    return visc;
}

void Vertex::printVertex()
{
    cout << "myid=" << id << " "; 
    cout << "nneigh= " << nneigh << " ntrian=" << ntrian << " ";
    cout << " x=" << xyz.x << " y=" << xyz.y << " z=" << xyz.z << " : ";
    cout << " x=" << force.x << " y=" << force.y << " z=" << force.z << " : ";
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