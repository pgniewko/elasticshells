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
        bonded_vertices.push_back( orig.bonded_vertices[i] );
    }

    for (int i = 0; i < facets_number; i++)
    {
        bonded_elements.push_back( orig.bonded_elements[i] );
    }
}

Vertex::~Vertex() {}

void Vertex::add_neighbor(int idx)
{
    try
    {
        if (idx < 0)
            throw RunTimeError("Trying to add a vertex with a negative index.\n"
                               "Runtime data is incorrect. Simulation will be terminated.\n");


        for (int i = 0; i < vertex_degree; i++)
        {
            if (bonded_vertices[i] == idx)
            {
                return;
            }
        }

        bonded_vertices.push_back( idx );
        vertex_degree++;
    }
    catch (RunTimeError& e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Vertex::add_element(int idx)
{
    try
    {
        if (idx < 0)
            throw RunTimeError("Trying to add a triangle with a negative index.\n"
                               "Runtime data is incorrect. Simulation will be terminated.\n");

        for (int i = 0; i < facets_number; i++)
        {
            if (bonded_elements[i] == idx)
            {
                return;
            }
        }

        bonded_elements.push_back( idx );
        facets_number++;
    }
    catch (RunTimeError& e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}

bool Vertex::is_neighbor(int vidx) const
{
    for (int i = 0; i < vertex_degree; i++)
    {
        if (bonded_vertices[i] == vidx)
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
    return bonded_vertices[idx];
}

int Vertex::getTriangleId(int idx) const
{
    return bonded_elements[idx];
}

std::ostream& operator<< (std::ostream& out, const Vertex& v)
{
    out << ' ' << v.my_id << ' ' << v.vertex_degree << ' ' << v.facets_number << ' ';

    for (int i = 0; i < v.vertex_degree; i++)
    {
        out << v.bonded_vertices[i] << ' ';
    }

    for (int i = 0; i < v.facets_number; i++)
    {
        out << v.bonded_elements[i] << ' ';
    }

    return out;
}