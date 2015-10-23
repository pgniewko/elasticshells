#include "Vertex.h"

Vertex::Vertex() : r_c(0, 0, 0), numBonded(0), numTris(0), numNbNeighs(0),
    domainIdx(-1), myid(-1), mass(1.0), visc(100.0), gtimer(0.0), my_type(vertex_t::MOTHER) {}

Vertex::Vertex(double x, double y, double z) : r_c(x, y, z), numBonded(0), numTris(0),
    numNbNeighs(0), domainIdx(-1), myid(-1), mass(1.0), visc(100.0), gtimer(0.0), my_type(vertex_t::MOTHER) {}

Vertex::Vertex(Vector3D v) : r_c(v), numBonded(0), numTris(0),
    numNbNeighs(0), domainIdx(-1), myid(-1), mass(1.0), visc(100.0), gtimer(0.0), my_type(vertex_t::MOTHER) {}

Vertex::Vertex(const Vertex& orig) : r_c(orig.r_c), f_c(orig.f_c), v_c(orig.v_c),
    r_p(orig.r_p), v_p(orig.v_p), f_p(orig.f_p),
    numBonded(orig.numBonded), numTris(orig.numTris), numNbNeighs(orig.numNbNeighs), domainIdx(orig.domainIdx),
    myid(orig.myid), mass( orig.mass ), visc(orig.visc), gtimer(0.0), my_type(orig.my_type)
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
        if (numBonded >= NEIGH_MAX)
            throw MaxSizeException("Maximum number of neighbors has been reached."
                                   "New neighbor will not be added!\n"
                                   "Simulation will be terminated.\n");

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
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (RunTimeError& e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Vertex::removeNeighbor(int vidx)
{
    int pos = -1;

    if ( isNeighbor(vidx) )
    {
        for (int i = 0; i < numBonded; i++)
        {
            if (bondedVerts[i] == vidx)
            {
                pos = i;
            }
        }

        for (int i = pos; i < numBonded - 1; i++)
        {
            bondedVerts[i] = bondedVerts[i + 1];
            r0[i] = r0[i + 1];
        }

        numBonded--;
    }

    return;
}

void Vertex::addTriangle(int idx)
{
    try
    {
        if (numTris >= TRIAN_MAX)
            throw MaxSizeException("Maximum number of triangles has been reached."
                                   "New triangle will not be added!\n"
                                   "Simulation will be terminated.\n");

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
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (RunTimeError& e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
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
    {
        return;
    }

    for (int i = pos; i < numTris - 1; i++)
    {
        bondedTris[i] = bondedTris[i + 1];
    }

    numTris--;
    return;
}

void Vertex::addNbNeighbor(int vertIdx, int cellIdx)
{
    try
    {
        if (numNbNeighs >= NBNEI_MAX)
            throw MaxSizeException("Maximum number of neighbors has been reached."
                                   "New neighbor will not be added.\n"
                                   "Simulation will be terminated!\n");

        if (vertIdx < 0)
            throw RunTimeError("Trying to add a vertex with a negative index.\n"
                               "Runtime data is incorrect. Simulation will be terminated!\n");

        if (cellIdx < 0)
            throw RunTimeError("Trying to add a vertex with a negative cell index.\n"
                               "Runtime data is incorrect. Simulation will be terminated!\n");

        nbVerts[numNbNeighs] = vertIdx;
        nbCellsIdx[numNbNeighs] = cellIdx;
        numNbNeighs++;
    }
    catch (MaxSizeException& e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (RunTimeError& e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}


bool Vertex::isNeighbor(int vidx) const
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
    f_c.x = 0.0;
    f_c.y = 0.0;
    f_c.z = 0.0;
}

void Vertex::voidVelocity()
{
    v_c.x = 0.0;
    v_c.y = 0.0;
    v_c.z = 0.0;
}

int Vertex::setId(int idx)
{
    myid = idx;
    return myid;
}

int Vertex::getId() const
{
    return myid;
}

int Vertex::getNumNeighs() const
{
    return numBonded;
}

int Vertex::getNumTris() const
{
    return numTris;
}

int Vertex::getNeighborId(int idx) const
{
    return bondedVerts[idx];
}

int Vertex::getTriangleId(int idx) const
{
    return bondedTris[idx];
}

double Vertex::getNeighborR0(int idx) const
{
    return r0[idx];
}

double Vertex::setMass(double m)
{
    mass = m;
    return mass;
}

double Vertex::getMass() const
{
    return mass;
}

double Vertex::setVisc(double v)
{
    visc = v;
    return visc;
}

double Vertex::getVisc() const
{
    return visc;
}

void Vertex::printVertex()
{
    std::cout << "myid=" << myid << " ";
    std::cout << "nneigh= " << numBonded << " ntrian=" << numTris << " ";
    std::cout << " x=" << r_c.x << " y=" << r_c.y << " z=" << r_c.z << " : ";
    std::cout << " x=" << f_c.x << " y=" << f_c.y << " z=" << f_c.z << " : ";

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

const vertex_t& Vertex::getMyType() const
{
    return my_type;
}

void Vertex::addTime(double dt)
{
    gtimer += dt;
}

void Vertex::voidTime()
{
    gtimer = 0.0;
}