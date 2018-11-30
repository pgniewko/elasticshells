#include "Vertex.h"

Vertex::Vertex() : r_c(0, 0, 0), r_p(0, 0, 0), vertex_degree(0), facets_number(0), myid(-1), myCellId(-1)
{
}

Vertex::Vertex(double x, double y, double z) : r_c(x, y, z), r_p(x, y, z), vertex_degree(0), facets_number(0),
    myid(-1), myCellId(-1)
{
}

Vertex::Vertex(const Vertex& orig) : r_c(orig.r_c), f_c(orig.f_c), r_p(orig.r_p), f_p(orig.f_p),
    v_c(orig.v_c),
    vertex_degree(orig.vertex_degree), facets_number(orig.facets_number), myid(orig.myid), myCellId(orig.myCellId)
{
    for (int i = 0; i < vertex_degree; i++)
    {
        bondedVerts.push_back( orig.bondedVerts[i] );
    }

    for (int i = 0; i < facets_number; i++)
    {
        bondedTris.push_back( orig.bondedTris[i] );
    }
}

Vertex::~Vertex() 
{
}

void Vertex::addNeighbor(int idx)
{
    try
    {
        if (idx < 0)
            throw RunTimeError("Trying to add a vertex with a negative index.\n"
                               "Runtime data is incorrect. Simulation will be terminated.\n");


        for (int i = 0; i < vertex_degree; i++)
        {
            if (bondedVerts[i] == idx)
            {
                return;
            }
        }

        bondedVerts.push_back( idx );
        vertex_degree++;
    }
    catch (RunTimeError& e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Vertex::addTriangle(int idx)
{
    try
    {
        if (idx < 0)
            throw RunTimeError("Trying to add a triangle with a negative index.\n"
                               "Runtime data is incorrect. Simulation will be terminated.\n");

        for (int i = 0; i < facets_number; i++)
        {
            if (bondedTris[i] == idx)
            {
                return;
            }
        }

        bondedTris.push_back( idx );
        facets_number++;
    }
    catch (RunTimeError& e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}

bool Vertex::isNeighbor(int vidx) const
{
    for (int i = 0; i < vertex_degree; i++)
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

int Vertex::set_id(int idx)
{
    myid = idx;
    return myid;
}

int Vertex::get_id() const
{
    return myid;
}

int Vertex::set_shell_id(int cellId)
{
    myCellId = cellId;
    return myCellId;
}

int Vertex::get_shell_id() const
{
    return myCellId;
}

int Vertex::get_vertex_degree() const
{
    return vertex_degree;
}

int Vertex::getNumTris() const
{
    return facets_number;
}

int Vertex::getNeighborId(int idx) const
{
    return bondedVerts[idx];
}

int Vertex::getTriangleId(int idx) const
{
    return bondedTris[idx];
}

// TODO:
// zamien te funkcje na przeciazony >> operator
// http://www.learncpp.com/cpp-tutorial/93-overloading-the-io-operators/

void Vertex::printVertex() const
{
    std::cout << "myid=" << myid << " ";
    std::cout << "nneigh= " << vertex_degree << " ntrian=" << facets_number << " ";
    std::cout << " x=" << r_c.x << " y=" << r_c.y << " z=" << r_c.z << " : ";
    std::cout << " x=" << f_c.x << " y=" << f_c.y << " z=" << f_c.z << " : ";

    for (int i = 0; i < vertex_degree; i++)
    {
        std::cout << bondedVerts[i] << " " ;
    }

    std::cout << " : ";

    for (int i = 0; i < facets_number; i++)
    {
        std::cout << bondedTris[i] << " " ;
    }

    std::cout << std::endl;
    return;
}

std::ostream& operator<< (std::ostream& out, const Vertex& v)
{
    out << ' ' << v.myid << ' ' << v.vertex_degree << ' ' << v.facets_number << ' ';

    for (int i = 0; i < v.vertex_degree; i++)
    {
        out << v.bondedVerts[i] << ' ';
    }

    for (int i = 0; i < v.facets_number; i++)
    {
        out << v.bondedTris[i] << ' ';
    }

    return out;
}