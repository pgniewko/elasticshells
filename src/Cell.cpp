#include "Cell.h"

//Cell::Cell() : numberV(0), numberT(0)
//{
//}

Cell::Cell(list<Triangle> tris) : numberV(0), numberT(0)
{
    constructVertices(tris);
    constructVTriangles(tris);
}

Cell::Cell(const Cell& orig) {
}

Cell::~Cell() {
}

void Cell::constructVertices(list<Triangle> tris)
{
    
    list<Vector3D> vectors;
    double xtmp, ytmp, ztmp;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
    {
        if ( isUnique(vectors, i->a) )
        {
            vectors.push_back(i->a);
            xtmp = i->a.x;
            ytmp = i->a.y;
            ztmp = i->a.z;
            vertices[numberV] = Vertex(xtmp, ytmp, ztmp);
            vertices[numberV].setId(numberV);
            numberV++;
        }
        
        if ( isUnique(vectors, i->b) )
        {
            vectors.push_back(i->b);
            xtmp = i->b.x;
            ytmp = i->b.y;
            ztmp = i->b.z;
            vertices[numberV] = Vertex(xtmp, ytmp, ztmp);
            vertices[numberV].setId(numberV);
            numberV++;
        }
        
        if ( isUnique(vectors, i->c) )
        {
            vectors.push_back(i->c);
            xtmp = i->c.x;
            ytmp = i->c.y;
            ztmp = i->c.z;
            vertices[numberV] = Vertex(xtmp, ytmp, ztmp);
            vertices[numberV].setId(numberV);
            numberV++;
        }
        
    }
}

void Cell::constructVTriangles(list<Triangle> tris)
{
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
    {
        Vertex *ptra;
        Vertex *ptrb;
        Vertex *ptrc;
        ptra = getVertex(i->a);
        ptrb = getVertex(i->b);
        ptrc = getVertex(i->c);
        VertexTriangle vrxt(ptra, ptrb, ptrc);
        triangles[numberT] = VertexTriangle(vrxt);
        numberT++;
    }
}

Vertex * Cell::getVertex(const Vector3D v)
{
    for (int i = 0; i < numberV; i++)
    {
        if (vertices[i].xyz.x == v.x && vertices[i].xyz.y == v.y && vertices[i].xyz.z == v.z)
        {
            return &vertices[i];
        }
    }
    return NULL;
 }


bool Cell::isUnique(list<Vector3D>& vlist, Vector3D& v)
{
    for(std::list<Vector3D>::iterator i = vlist.begin(); i != vlist.end(); ++i) 
    {
        if (i->x == v.x && i->y == v.y && i->z == v.z)
        {
            return false;
        }
    }
    return true;
}

double Cell::calcSurfaceArea()
{
    double surface = 0.0;
    for (int i = 0; i < numberT; i++)
    {
        surface += triangles[i].area();
    }

    return surface;
}

double Cell::calcVolume()
{
    double volume = 0.0;
    for (int i = 0; i < numberT; i++)
    {
        Tetrahedron tetr(triangles[i].a->xyz, triangles[i].b->xyz, triangles[i].c->xyz, cm);
        volume += tetr.volume();
    }

    return volume;
}

void Cell::calcCM()
{
    Vector3D tmp(0,0,0);
    
    double M = 0.0;
    for (int i = 0; i < numberV; i++)
    {
        tmp += vertices[i].xyz;
        M += vertices[i].getMass();
    }

    tmp /= M;
    cm = tmp;
}

int Cell::numberofFaces() 
{
    return numberT;
}

int Cell::numberofVertices() 
{
    return numberV;
}