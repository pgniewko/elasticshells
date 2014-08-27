#ifndef VERTEXTRIANGLE_H
#define	VERTEXTRIANGLE_H

#include "Vertex.h"
#include "Triangle.h"

class Vertex; 

class VertexTriangle {
public:
    VertexTriangle();
    VertexTriangle(int, int, int);
    VertexTriangle(const VertexTriangle& orig);
    virtual ~VertexTriangle();
    void setId(int);
    double area(const Vertex[]);
    
    void printVertexTriangle();
    
    int ia, ib, ic;
    int id;
};

#endif	/* VERTEXTRIANGLE_H */
