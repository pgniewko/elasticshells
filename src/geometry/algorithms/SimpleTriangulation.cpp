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

list<Triangle> SimpleTriangulation::triangulate()
{
    if (depth < 0)
    {
        return NULL;
    }
    
    else if (depth == 0)
    {
       tris.push_back( Triangle(Point(0,0,0), Point(0,0,0), Point(0,0,0) ) );
    } 
    else if (depth == 1)
    {
        createCube();
        for (int i = 0; i < (depth-1); i++)
        {
            subdivide();
        }
    }
    
    return tris;
}

void SimpleTriangulation::createCube()
{
    // create cube
    //top
    tris.push_back(Triangle(Point(1, 1, 1), Point(-1, 1, -1), Point(1, 1, -1)));
    tris.push_back(Triangle(Point(-1, 1, -1), Point(1, 1, 1), Point(-1, 1, 1)));
    //bottom
    tris.push_back(Triangle(Point(1, -1, 1), Point(-1, -1, -1), Point(1, -1, -1)));
    tris.push_back(Triangle(Point(-1, -1, -1), Point(1, -1, 1), Point(-1, -1, 1)));

    //right
    tris.push_back(Triangle(Point(1, 1, 1), Point(1, -1, -1), Point(1, 1, -1)));
    tris.push_back(Triangle(Point(1, -1, -1), Point(1, 1, 1), Point(1, -1, 1)));
    //left
    tris.push_back(Triangle(Point(-1, 1, 1), Point(-1, -1, -1), Point(-1, 1, -1)));
    tris.push_back(Triangle(Point(-1, -1, -1), Point(-1, 1, 1), Point(-1, -1, 1)));

    //front
    tris.push_back(Triangle(Point(1, 1, 1), Point(-1, -1, 1), Point(1, -1, 1)));
    tris.push_back(Triangle(Point(-1, -1, 1), Point(1, 1, 1), Point(-1, 1, 1)));
    //rear
    tris.push_back(Triangle(Point(1, 1, -1), Point(-1, -1, -1), Point(1, -1, -1)));
    tris.push_back(Triangle(Point(-1, -1, -1), Point(1, 1, -1), Point(-1, 1, -1)));
    /**/
}

void SimpleTriangulation::subdivide()
{
//    printf("jestem ... \n");  
    int counter = 0;
    std::list<Triangle> newTris;
    float l = tris.begin()->a.length();
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) { // go through all triangles
        Point mid = (i->a + i->b) * 0.5f; // point betwenn points A and B
        mid.setLength(l); // put in on the sphere
        newTris.push_back(Triangle(i->b, i->c, mid)); // remember new triangles
        newTris.push_back(Triangle(i->a, i->c, mid));
//        counter++;
//        printf("trojkat no: %d\n", counter);
    }
    tris.swap(newTris); // use new set of triangles;
}