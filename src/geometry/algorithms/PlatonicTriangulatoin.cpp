#include "PlatonicTriangulatoin.h"

PlatonicTriangulatoin::PlatonicTriangulatoin() : depth(0), type(0)
{}

PlatonicTriangulatoin::PlatonicTriangulatoin(int d) : depth(0), type(0)
{
    depth = (d >= 0) ? d : 0;
}

PlatonicTriangulatoin::PlatonicTriangulatoin(int d, int t) : depth(0), type(0)
{
    depth = (d >= 0) ? d : 0;
    type  = (t >= 0) ? t : 0;
}

PlatonicTriangulatoin::PlatonicTriangulatoin(const PlatonicTriangulatoin& orig) :
    depth(orig.depth), type(orig.type) {}

PlatonicTriangulatoin::~PlatonicTriangulatoin() {}

std::list<Triangle> PlatonicTriangulatoin::triangulate(double r0)
{
    if (depth == 0)
    {
        tris.push_back( Triangle(Vector3D(0, 0, 0), Vector3D(0, 0, 0), Vector3D(0, 0, 0)) );
    }
    else
    {
        if (type == 0)
        {
            createTetrahedron();
        }
        else if (type == 1)
        {
            createHexahedron();
        }
        else if (type == 2)
        {
            createOctahedron();
        }
        else if (type == 3)
        {
            createIcosahedron();
        }
        else
        {
            createTetrahedron();
        }

        for (int i = 0; i < (depth - 1); i++)
        {
            subdivide();
        }
    }

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        i->a.set_length(r0);
        i->b.set_length(r0);
        i->c.set_length(r0);
    }

    return tris;
}


std::list<Triangle> PlatonicTriangulatoin::triangulate()
{
    return triangulate(1.0);
}


void PlatonicTriangulatoin::createTetrahedron()
{
    double C = constants::sqrt3 / 3.0; // - OK, checked. change C = SQRT3/3 to 1. no diff.
    tris.push_back(Triangle(Vector3D( C, C, C), Vector3D(-C, -C, C), Vector3D(-C, C, -C) ));
    tris.push_back(Triangle(Vector3D( C, C, C), Vector3D( C, -C, -C), Vector3D(-C, -C, C) ));
    tris.push_back(Triangle(Vector3D(-C, C, -C), Vector3D(-C, -C, C), Vector3D( C, -C, -C) ));
    tris.push_back(Triangle(Vector3D( C, -C, -C), Vector3D( C, C, C), Vector3D(-C, C, -C) ));
}

void PlatonicTriangulatoin::createHexahedron()
{
    // create cube - OK, checked
    //top
    tris.push_back(Triangle(Vector3D( 1, 1, 1), Vector3D(-1, 1, -1), Vector3D( 1, 1, -1)));
    tris.push_back(Triangle(Vector3D(-1, 1, -1), Vector3D( 1, 1, 1), Vector3D(-1, 1, 1)));
    //bottom
    tris.push_back(Triangle(Vector3D( 1, -1, 1), Vector3D(-1, -1, -1), Vector3D( 1, -1, -1)));
    tris.push_back(Triangle(Vector3D(-1, -1, -1), Vector3D( 1, -1, 1), Vector3D(-1, -1, 1)));
    //right
    tris.push_back(Triangle(Vector3D( 1, 1, 1), Vector3D( 1, -1, -1), Vector3D( 1, 1, -1)));
    tris.push_back(Triangle(Vector3D( 1, -1, -1), Vector3D( 1, 1, 1), Vector3D( 1, -1, 1)));
    //left
    tris.push_back(Triangle(Vector3D(-1, 1, 1), Vector3D(-1, -1, -1), Vector3D(-1, 1, -1)));
    tris.push_back(Triangle(Vector3D(-1, -1, -1), Vector3D(-1, 1, 1), Vector3D(-1, -1, 1)));
    //front
    tris.push_back(Triangle(Vector3D( 1, 1, 1), Vector3D(-1, -1, 1), Vector3D( 1, -1, 1)));
    tris.push_back(Triangle(Vector3D(-1, -1, 1), Vector3D( 1, 1, 1), Vector3D(-1, 1, 1)));
    //rear
    tris.push_back(Triangle(Vector3D( 1, 1, -1), Vector3D(-1, -1, -1), Vector3D( 1, -1, -1)));
    tris.push_back(Triangle(Vector3D(-1, -1, -1), Vector3D( 1, 1, -1), Vector3D(-1, 1, -1)));
    /**/
}

void PlatonicTriangulatoin::createOctahedron()
{
    // Join vertices to create a unit octahedron - OK, checked
    // The top half
    tris.push_back(Triangle(Vector3D( 1, 0, 0), Vector3D( 0, 0, 1), Vector3D( 0, 1, 0)));
    tris.push_back(Triangle(Vector3D( 0, 1, 0), Vector3D( 0, 0, 1), Vector3D(-1, 0, 0)));
    tris.push_back(Triangle(Vector3D(-1, 0, 0), Vector3D( 0, 0, 1), Vector3D( 0, -1, 0)));
    tris.push_back(Triangle(Vector3D( 0, -1, 0), Vector3D( 0, 0, 1), Vector3D( 1, 0, 0)));
    // The bottom half
    tris.push_back(Triangle(Vector3D( 1, 0, 0), Vector3D( 0, 1, 0), Vector3D( 0, 0, -1)));
    tris.push_back(Triangle(Vector3D( 0, 1, 0), Vector3D(-1, 0, 0), Vector3D( 0, 0, -1)));
    tris.push_back(Triangle(Vector3D(-1, 0, 0), Vector3D( 0, -1, 0), Vector3D( 0, 0, -1)));
    tris.push_back(Triangle(Vector3D( 0, -1, 0), Vector3D( 1, 0, 0), Vector3D( 0, 0, -1)));
}

void PlatonicTriangulatoin::createIcosahedron()
{
    //Twelve vertices of icosahedron on unit sphere
    double t = 0.8506508084; // t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)
    double u = 0.5257311121; // one=1/sqrt(1+t^2) , unit sphere
    // Structure for unit icosahedron - OK, checked
    tris.push_back(Triangle(Vector3D( u, 0, t), Vector3D(-u, 0, t), Vector3D( 0, t, u)));
    tris.push_back(Triangle(Vector3D( u, 0, t), Vector3D( 0, -t, u), Vector3D(-u, 0, t)));
    tris.push_back(Triangle(Vector3D( u, 0, -t), Vector3D( 0, t, -u), Vector3D(-u, 0, -t)));
    tris.push_back(Triangle(Vector3D( u, 0, -t), Vector3D(-u, 0, -t), Vector3D( 0, -t, -u)));
    tris.push_back(Triangle(Vector3D( t, u, 0), Vector3D( t, -u, 0), Vector3D( u, 0, t)));
    tris.push_back(Triangle(Vector3D( t, u, 0), Vector3D( u, 0, -t), Vector3D( t, -u, 0)));
    tris.push_back(Triangle(Vector3D(-t, -u, 0), Vector3D(-t, u, 0), Vector3D(-u, 0, t)));
    tris.push_back(Triangle(Vector3D(-t, -u, 0), Vector3D(-u, 0, -t), Vector3D(-t, u, 0)));
    tris.push_back(Triangle(Vector3D( 0, t, u), Vector3D( 0, t, -u), Vector3D( t, u, 0)));
    tris.push_back(Triangle(Vector3D( 0, t, u), Vector3D(-t, u, 0), Vector3D( 0, t, -u)));
    tris.push_back(Triangle(Vector3D( 0, -t, u), Vector3D( t, -u, 0), Vector3D( 0, -t, -u)));
    tris.push_back(Triangle(Vector3D( 0, -t, u), Vector3D( 0, -t, -u), Vector3D(-t, -u, 0)));
    tris.push_back(Triangle(Vector3D( 0, t, u), Vector3D( t, u, 0), Vector3D( u, 0, t)));
    tris.push_back(Triangle(Vector3D( 0, t, -u), Vector3D( u, 0, -t), Vector3D( t, u, 0)));
    tris.push_back(Triangle(Vector3D( u, 0, t), Vector3D( t, -u, 0), Vector3D( 0, -t, u)));
    tris.push_back(Triangle(Vector3D( u, 0, -t), Vector3D( 0, -t, -u), Vector3D( t, -u, 0)));
    tris.push_back(Triangle(Vector3D(-u, 0, t), Vector3D(-t, u, 0), Vector3D( 0, t, u)));
    tris.push_back(Triangle(Vector3D(-u, 0, -t), Vector3D( 0, t, -u), Vector3D(-t, u, 0)));
    tris.push_back(Triangle(Vector3D(-u, 0, t), Vector3D( 0, -t, u), Vector3D(-t, -u, 0)));
    tris.push_back(Triangle(Vector3D(-u, 0, -t), Vector3D(-t, -u, 0), Vector3D( 0, -t, -u)));
}

/*
 * Subdivide each triangle in the old approximation and normalize
 * the new points thus generated to lie on the surface of the unit
 * sphere.
 * Each input triangle with vertices labeled [0,1,2] as shown
 * below will be turned into four new triangles:
 *
 *             Make new points
 *                  a = (0+2)/2
 *                  b = (0+1)/2
 *                  c = (1+2)/2
 *         1
 *        /\        Normalize a, b, c
 *       /  \
 *     b/____\ c    Construct new triangles
 *     /\    /\       t1 [0,b,a]
 *    /  \  /  \      t2 [b,1,c]
 *   /____\/____\     t3 [a,b,c]
 *  0      a     2    t4 [a,c,2]
 *
 */
void PlatonicTriangulatoin::subdivide()
{
    int counter = 0;
    std::list<Triangle> newTris;
    float l = tris.begin()->a.length();

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)  // go through all triangles
    {
        Vector3D abh = (i->a + i->b) * 0.5; // point between points A and B
        Vector3D ach = (i->a + i->c) * 0.5; // point between points A and C
        Vector3D bch = (i->b + i->c) * 0.5; // point between points B and C
        abh.set_length(l);
        ach.set_length(l);
        bch.set_length(l);
        newTris.push_back(Triangle(i->a, abh, ach));
        newTris.push_back(Triangle(ach, i->c, bch));
        newTris.push_back(Triangle(abh, bch, i->b));
        newTris.push_back(Triangle(abh, ach, bch));
        counter++;
    }

    tris.swap(newTris); // use new set of triangles;
}