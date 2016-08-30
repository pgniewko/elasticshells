#include "Vertex.h"

Vertex::Vertex() : r_c(0, 0, 0), r_p(0, 0, 0), numBonded(0), numTris(0), myid(-1), myCellId(-1)
{
    for (int i = 0; i < NEIGH_MAX; i++)
    {
        bondedVerts[i] = -1;
        r0[i] = 0.0;
        k0[i] = 0.0;
    }

    for (int i = 0; i < TRIAN_MAX; i++)
    {
        bondedTris[i] = -1;
    }
}

Vertex::Vertex(double x, double y, double z) : r_c(x, y, z), r_p(x, y, z), numBonded(0), numTris(0),
    myid(-1), myCellId(-1)
{
    for (int i = 0; i < NEIGH_MAX; i++)
    {
        bondedVerts[i] = -1;
        r0[i] = 0.0;
        k0[i] = 0.0;
    }

    for (int i = 0; i < TRIAN_MAX; i++)
    {
        bondedTris[i] = -1;
    }
}

Vertex::Vertex(Vector3D v) : r_c(v), r_p(v), numBonded(0), numTris(0), myid(-1), myCellId(-1)
{
    for (int i = 0; i < NEIGH_MAX; i++)
    {
        bondedVerts[i] = -1;
        r0[i] = 0.0;
        k0[i] = 0.0;
    }

    for (int i = 0; i < TRIAN_MAX; i++)
    {
        bondedTris[i] = -1;
    }
}

Vertex::Vertex(const Vertex& orig) : r_c(orig.r_c), f_c(orig.f_c), r_p(orig.r_p), f_p(orig.f_p),
    v_p(orig.v_p), v_c(orig.v_c), a_p(orig.a_p), a_c(orig.a_c),
    numBonded(orig.numBonded), numTris(orig.numTris), myid(orig.myid), myCellId(orig.myCellId)
{
    for (int i = 0; i < numBonded; i++)
    {
        bondedVerts[i] = orig.bondedVerts[i];
        r0[i] = orig.r0[i];
        k0[i] = orig.k0[i];
    }

    for (int i = 0; i < numTris; i++)
    {
        bondedTris[i] = orig.bondedTris[i];
    }
}

Vertex::~Vertex() {}

void Vertex::addNeighbor(int idx, double r0v)
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
        r0[numBonded] = r0v;
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

void Vertex::voidForce()
{
    f_c.x = 0.0;
    f_c.y = 0.0;
    f_c.z = 0.0;
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

int Vertex::setCellId(int cellId)
{
    myCellId = cellId;
    return myCellId;
}

int Vertex::getCellId() const
{
    return myCellId;
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

// TODO:
// zamien te funkcje na przeciazony >> operator
// http://www.learncpp.com/cpp-tutorial/93-overloading-the-io-operators/

void Vertex::printVertex() const
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
    return;
}

std::ostream& operator<< (std::ostream& out, const Vertex& v)
{
    out << ' ' << v.myid << ' ' << v.numBonded << ' ' << v.numTris << ' ';

    for (int i = 0; i < v.numBonded; i++)
    {
        out << v.bondedVerts[i] << ' ' << v.r0[i] << ' ' << v.k0[i] << ' ';
    }

    for (int i = 0; i < v.numTris; i++)
    {
        out << v.bondedTris[i] << ' ';
    }

    return out;
}