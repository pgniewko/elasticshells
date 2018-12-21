#include "SimpleTriangulation.h"

SimpleTriangulation::SimpleTriangulation()
{
    depth = 0;
}

SimpleTriangulation::SimpleTriangulation(int d)
{
    depth = (d >= 0) ? d : 0;
}

SimpleTriangulation::SimpleTriangulation(const SimpleTriangulation& orig) : depth(orig.depth) {}

SimpleTriangulation::~SimpleTriangulation() {}

std::list<Triangle> SimpleTriangulation::triangulate(double r0)
{
//    if (depth < 0)
//    {
//        return NULL;
    if (depth == 0)
    {
        tris.push_back( Triangle(Vector3D(0, 0, 0), Vector3D(0, 0, 0), Vector3D(0, 0, 0)) );
    }
    else
    {
        create_cube();

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

std::list<Triangle> SimpleTriangulation::triangulate()
{
    return triangulate(1.0);
}

void SimpleTriangulation::create_cube()
{
    // create cube
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

void SimpleTriangulation::subdivide()
{
    int counter = 0;
    std::list<Triangle> newTris;
    float l = tris.begin()->a.length();

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)  // go through all triangles
    {
        Vector3D mid = (i->a + i->b) * 0.5f; // point between points A and B
        mid.set_length(l); // put in on the sphere
        newTris.push_back(Triangle(i->b, i->c, mid)); // remember new triangles
        newTris.push_back(Triangle(i->a, i->c, mid));
        counter++;
    }

    tris.swap(newTris); // use new set of triangles;
}