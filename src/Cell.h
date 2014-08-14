#ifndef CELL_H
#define	CELL_H

#define MAX_V 200
#define MAX_T 500

#include <list>

#include "geometry/Vector3D.h"
#include "geometry/Tetrahedron.h"
#include "geometry/Triangle.h"
#include "geometry/Vertex.h"
#include "geometry/VertexTriangle.h"

#include "force/HookeanForce.h"
#include "force/OsmoticForce.h"
#include "force/NbRepulsiveForce.h"

using namespace std;

class Cell {
public:
//    Cell();
    Cell(list<Triangle>);
    Cell(const Cell& orig);
    virtual ~Cell();
    double calcSurfaceArea();
    double calcVolume();
    void calcCM();
    int numberofFaces() ;
    int numberofVertices();
    
    
//    void createDataStructure();
    void setTopology();
    void printTopology();
    void calcForces();
    void calcForces(const Cell&);
    void addVelocity(const Vector3D&);
    void addXYZ(const Vector3D&);
    Vector3D cm;
    Vertex vertices[MAX_V];
    VertexTriangle triangles[MAX_T];
    int numberV;
    int numberT;
    
    void saveTriangulatedSurface(const char*);
    void saveRenderingScript(const char*, const char*);
    
    
private:
    bool isUnique(list<Vector3D>&, Vector3D&);
    Vertex * getVertex(const Vector3D);
    void constructVertices(list<Triangle>);
    void constructVTriangles(list<Triangle>);
    

};

#endif	/* CELL_H */