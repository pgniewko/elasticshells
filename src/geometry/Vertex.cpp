#include "Vertex.h"

Vertex::Vertex() : r_c(0.0, 0.0, 0.0), 
        vertex_degree(0), 
        facets_number(0), 
        my_id(-1), 
        my_shell_id(-1)
{
}

Vertex::Vertex(double x, double y, double z) : r_c(x, y, z), 
        vertex_degree(0), 
        facets_number(0),
        my_id(-1), my_shell_id(-1)
{
}

Vertex::Vertex(const Vertex& orig) : r_c(orig.r_c),
        vertex_degree(orig.vertex_degree), 
        facets_number(orig.facets_number), 
        my_id(orig.my_id), 
        my_shell_id(orig.my_shell_id)
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

Vertex::~Vertex() {}

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

int Vertex::set_id(int idx)
{
    my_id = idx;
    return my_id;
}

int Vertex::get_id() const
{
    return my_id;
}

int Vertex::set_shell_id(int cellId)
{
    my_shell_id = cellId;
    return my_shell_id;
}

int Vertex::get_shell_id() const
{
    return my_shell_id;
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

std::ostream& operator<< (std::ostream& out, const Vertex& v)
{
    out << ' ' << v.my_id << ' ' << v.vertex_degree << ' ' << v.facets_number << ' ';

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