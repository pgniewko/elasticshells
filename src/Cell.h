#include <list>
#include "geometry/Triangle.h"

#ifndef CELL_H
#define	CELL_H

#include "geometry/Vector3D.h"
#include "geometry/Tetrahedron.h"
#include "geometry/Vertex.h"

using namespace std;

class Cell {
public:
//    Cell();
    Cell(list<Triangle>);
    Cell(const Cell& orig);
    virtual ~Cell();
    double surfaceArea();
    double volume();
    void calcCM();
    int numberofFaces() ;
    int numberofVertices();
//    void createDataStructure();
//    void calcForces();
//    void calcForces(const Cell&);
    Vector3D cm;
    list<Triangle> tris;
    list<Vertex> vertices;
private:
    bool isUnique(list<Vector3D>&, Vector3D&); 
//    void constructTriangles();
    

};

#endif	/* CELL_H */

