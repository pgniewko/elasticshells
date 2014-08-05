#include <list>
#include "geometry/Triangle.h"

#ifndef CELL_H
#define	CELL_H

#include "geometry/Vector3D.h"
#include "geometry/Tetrahedron.h"

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
    int numberofVertices() ;
    Vector3D cm;
    list<Triangle> tris;
private:
    bool isUnique(list<Vector3D>, Vector3D) ;
    

};

#endif	/* CELL_H */

