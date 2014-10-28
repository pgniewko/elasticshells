#include "Vertex.h"

Vertex::Vertex() : xyz(0, 0, 0), numBonded(0), numTris(0), numNbNeighs(0), 
        domainIdx(-1), myid(-1), mass(1.0) {}

Vertex::Vertex(double x, double y, double z) : xyz(x, y, z), numBonded(0), numTris(0), 
        numNbNeighs(0), domainIdx(-1), myid(-1), mass(1.0) {}

Vertex::Vertex(Vector3D v) : xyz(v), numBonded(0), numTris(0), 
        numNbNeighs(0), domainIdx(-1), myid(-1), mass(1.0) {}

Vertex::Vertex(const Vertex& orig) : xyz(orig.xyz), force(orig.force), velocity(orig.velocity),
    tmp_xyz(orig.tmp_xyz), tmp_force(orig.tmp_force), tmp_velocity(orig.tmp_velocity),
    numBonded(orig.numBonded), numTris(orig.numTris), numNbNeighs(orig.numNbNeighs), domainIdx(orig.domainIdx),
    myid(orig.myid), mass( orig.mass ), visc(orig.visc)
{
    for (int i = 0; i < numBonded; i++)
    {
        bondedVerts[i] = orig.bondedVerts[i];
        r0[i] = orig.r0[i];
    }

    for (int i = 0; i < numTris; i++)
    {
        bondedTris[i] = orig.bondedTris[i];
    }

    for (int i = 0; i < numNbNeighs; i++)
    {
        nbVerts[i] = orig.nbVerts[i];
        nbCellsIdx[i] = orig.nbCellsIdx[i];
    }
}

Vertex::~Vertex() {}

void Vertex::addNeighbor(int idx, double k0n)
{
    try
    {
        if (getNumNeighs() >= NEIGH_MAX)
            throw MaxSizeException("Maximum number of neighbors has been reached.\n"
                                   "New neighbor will not be added !\n");

        if (idx < 0)
            throw RunTimeError("Trying to add a vertex with a negative index.\n"
                    "Runtime data is incomplete. Simulation will be terminated. \n"); 
                
        for (int i = 0; i < numBonded; i++)
        {
            if (bondedVerts[i] == idx)
            {
                return;
            }
        }

        bondedVerts[numBonded] = idx;
        r0[numBonded] = k0n;
        numBonded++;
    }
    catch (MaxSizeException& e)
    {
        std::cout << e.what() << std::endl;
        return;
    }
    catch (RunTimeError& e)
    {
        std::cout << e.what() << std::endl;
        exit(1);
    }
}

void Vertex::addTriangle(int idx)
{
    try
    {
        if (getNumNeighs() >= TRIAN_MAX)
            throw MaxSizeException("Maximum number of triangles has been reached."
                                   "New triangle will not be added !");

        if (idx < 0)
            throw RunTimeError("Trying to add a triangle with a negative index.\n"
                    "Runtime data is incomplete. Simulation will be terminated. \n");        
        
        for (int i = 0; i < numTris; i++)
        {
            if (bondedTris[i] == idx)
            {
                return;
            }
        }

        bondedTris[numTris] = idx;
        numTris++;
    }
    catch (MaxSizeException& e)
    {
        std::cout << e.what() << std::endl;
        return;
    }
    catch (RunTimeError& e)
    {
        std::cout << e.what() << std::endl;
        exit(1);
    }
}

bool Vertex::isNeighbor(int vidx)
{
    for (int i = 0; i < numBonded; i++)
    {
        if (bondedVerts[i] == vidx)
        {
            return true;
        }
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
    myid = idx;
    return myid;
}

int Vertex::getId()
{
    return myid;
}

int Vertex::getNumNeighs()
{
    return numBonded;
}

int Vertex::getNumTris()
{
    return numTris;
}

int Vertex::getNeighborId(int idx)
{
    return bondedVerts[idx];
}

int Vertex::getTriangleId(int idx)
{
    return bondedTris[idx];
}

double Vertex::getNeighborR0(int idx)
{
    return r0[idx];
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
    std::cout << "myid=" << myid << " ";
    std::cout << "nneigh= " << numBonded << " ntrian=" << numTris << " ";
    std::cout << " x=" << xyz.x << " y=" << xyz.y << " z=" << xyz.z << " : ";
    std::cout << " x=" << force.x << " y=" << force.y << " z=" << force.z << " : ";

    for (int i = 0; i < numBonded; i++)
    {
        std::cout << bondedVerts[i] << " " ;
    }

    std::cout << " : ";

    for (int i = 0; i < numTris; i++)
    {
        std::cout << bondedTris[i] << " " ;
    }

    std::cout << std::endl;
}