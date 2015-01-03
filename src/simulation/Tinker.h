#ifndef TINKER_H
#define	TINKER_H

#include <list>
#include <vector>

#include "Environment.h"
#include "Cell.h"
#include "geometry/Triangle.h"
#include "geometry/Vertex.h"
#include "geometry/VertexTriangle.h"

class Cell;

class Tinker {
public:
    Tinker();
    Tinker(const Tinker& orig);
    virtual ~Tinker();
    static void constructVertices(Cell&, std::list<Triangle>&);
    static void constructVTriangles(Cell&, std::list<Triangle>&);
    static void constructTopology(Cell&);
    
    static void grow(Cell&);
    static void bud(Cell&);
    static void divide(Cell&);
private:
    static bool isUnique(std::list<Vector3D>&, Vector3D&);
    static int getRandomVertex(Cell&);
    static int getLeastDensiestVertex(Cell&);
    static int getVertex(Cell&, const Vector3D&);
    
    static int vidx[MAX_V];

};

#endif	/* TINKER_H */

