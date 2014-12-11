#include "Vertex.h"

Vertex::Vertex() : xyz(0, 0, 0), numBonded(0), numTris(0), numNbNeighs(0),
    domainIdx(-1), myid(-1), mass(1.0), isbud(false), my_type(vertex_t::MOTHER) {}

Vertex::Vertex(double x, double y, double z) : xyz(x, y, z), numBonded(0), numTris(0),
    numNbNeighs(0), domainIdx(-1), myid(-1), mass(1.0), isbud(false), my_type(vertex_t::MOTHER) {}

Vertex::Vertex(Vector3D v) : xyz(v), numBonded(0), numTris(0),
    numNbNeighs(0), domainIdx(-1), myid(-1), mass(1.0), isbud(false) {}

Vertex::Vertex(const Vertex& orig) : xyz(orig.xyz), force(orig.force), velocity(orig.velocity),
    tmp_xyz(orig.tmp_xyz), tmp_force(orig.tmp_force), tmp_velocity(orig.tmp_velocity),
    numBonded(orig.numBonded), numTris(orig.numTris), numNbNeighs(orig.numNbNeighs), domainIdx(orig.domainIdx),
    myid(orig.myid), mass( orig.mass ), visc(orig.visc), isbud(orig.isbud), my_type(orig.my_type)
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
    //std::cout <<" my index=" << myid << " adding idx=" << idx << std::endl;
    try
    {
        if (numBonded >= NEIGH_MAX)
            throw MaxSizeException("Maximum number of neighbors has been reached.\n"
                                   "New neighbor will not be added !\n");

        if (idx < 0)
            throw RunTimeError("Trying to add a vertex with a negative index.\n"
                               "Runtime data is incorrect. Simulation will be terminated.\n");

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

void Vertex::removeNeighbor(int vidx)
{
    //std::cout << "my index=" << myid << " removing vidx="<<vidx << std::endl;
    int pos;
    if ( isNeighbor(vidx) )
    {
        for (int i = 0; i < numBonded; i++)
        {
            if (bondedVerts[i] == vidx)
            {
                pos = i;
            }
        }
        
        //std::cout << "pos of vidx pos=" << pos << std::endl;
        
        for (int i = 0; i < numBonded ; i++)
        {
            //std::cout << " (bv, r0) = " << bondedVerts[i] << "," <<r0[i];
        }
        //std::cout << std::endl;
        
        for (int i = pos; i < numBonded - 1; i++)
        {
            bondedVerts[i] = bondedVerts[i + 1];
            r0[i] = r0[i + 1];
        }
        
        numBonded--;
        
        for (int i = 0; i < numBonded ; i++)
        {
            //std::cout << " (bv, r0) = " << bondedVerts[i] << "," <<r0[i];
        }
        //std::cout << std::endl;
    }
}

void Vertex::addTriangle(int idx)
{
    try
    {
        if (numTris >= TRIAN_MAX)
            throw MaxSizeException("Maximum number of triangles has been reached."
                                   "New triangle will not be added !");

        if (idx < 0)
            throw RunTimeError("Trying to add a triangle with a negative index.\n"
                               "Runtime data is incorrect. Simulation will be terminated.\n");

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

void Vertex::removeTriangle(int tidx)
{
    int pos = -1;
    for (int i = 0; i < numTris; i++)
    {
        if (bondedTris[i] == tidx)
        {
            pos = i;
        }
    }
    
    if (pos == -1)
        return;
        
    for (int i = pos; i < numTris - 1; i++)
    {
        bondedTris[i] = bondedTris[i + 1];
    }
    numTris--;
}

void Vertex::addNbNeighbor(int vertIdx, int cellIdx)
{
    try
    {
        if (numNbNeighs >= NBNEI_MAX)
            throw MaxSizeException("Maximum number of neighbors has been reached.\n"
                                   "New neighbor will not be added !\n");

        if (vertIdx < 0)
            throw RunTimeError("Trying to add a vertex with a negative index.\n"
                               "Runtime data is incorrect. Simulation will be terminated. \n");
        
        if (cellIdx < 0)
            throw RunTimeError("Trying to add a vertex with a negative cell index.\n"
                               "Runtime data is incorrect. Simulation will be terminated. \n");
        
        
        nbVerts[numNbNeighs] = vertIdx;
        nbCellsIdx[numNbNeighs] = cellIdx;
        numNbNeighs++;

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

void Vertex::sortNbList()
{
    std::vector<nblist_t> v_nb_lsit;
    
    for (int i = 0; i < numNbNeighs; i++)
    {
        struct nblist_t nblist;
        nblist.cell_id = nbCellsIdx[i];
        nblist.vertex_id = nbVerts[i];
	v_nb_lsit.push_back(nblist);
    }

    std::sort(v_nb_lsit.begin(), v_nb_lsit.end());
    
    for (unsigned int i = 0; i < v_nb_lsit.size(); i++)
    {
        nbCellsIdx[i] = v_nb_lsit[i].cell_id;
        nbVerts[i] = v_nb_lsit[i].vertex_id;
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

void Vertex::normalizedR0(double newR0)
{
    for (int i = 0; i < numBonded; i++)
    {
        r0[i] = newR0;
    }
}

bool Vertex::isBud()
{
    return isbud;
}

bool Vertex::setBud(bool boolval)
{
    isbud = boolval;
    return isbud;
}