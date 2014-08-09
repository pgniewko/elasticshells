#include "Cell.h"

//Cell::Cell() {
//}

Cell::Cell(list<Triangle> tlist) : tris(tlist) {}

Cell::Cell(const Cell& orig) {
}

Cell::~Cell() {
}

double Cell::surfaceArea()
{
    double surface = 0.0;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) {
        surface += i->area();
    }
    return surface;
}

double Cell::volume()
{
    double volume = 0.0;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) {
        Tetrahedron tetr(i->a, i->b, i->c, cm);
        volume += tetr.volume();
    }
    return volume;
}

void Cell::calcCM()
{
    Vector3D tmp(0,0,0);
    int counter = 0;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
    {
        tmp += i->a;
        tmp += i->b;
        tmp += i->b;
        counter += 3.0;
    }
    tmp /= counter;
    cm = tmp;
}

int Cell::numberofFaces() 
{
    return tris.size();
}

int Cell::numberofVertices() 
{
    
    list<Vector3D> vertices;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
    {
        cout << "next triangle" <<endl;
        if ( isUnique(vertices, i->a) )
            vertices.push_back(i->a);
        
        if ( isUnique(vertices, i->b) )
            vertices.push_back(i->b);
        
        if ( isUnique(vertices, i->c) )
            vertices.push_back(i->c);
         cout << "end of next triangle" <<endl;
    }
//    return 0;
    return vertices.size();
}

bool Cell::isUnique(list<Vector3D>& vlist, Vector3D& v)
{
    cout << "inside unique" << endl;
    for(std::list<Vector3D>::iterator i = vlist.begin(); i != vlist.end(); ++i) 
    {
        if (i->x == v.x && i->y == v.y && i->z == v.z)
        {
            return false;
        }
    }
    return true;
}


void Cell::createDataStructure()
{
    list<Vector3D> vertices;
    list<Triangle> triangles;
    int vertexId = 0;
    int triangleId = 0;
    
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
    {
        if ( isUnique(vertices, i->a) )
        {
            vertices.push_back(i->a);
        }
        
        if ( isUnique(vertices, i->b) )
            vertices.push_back(i->b);
        
        if ( isUnique(vertices, i->c) )
            vertices.push_back(i->c);
    }
//    return 0;
    return vertices.size();   
}

