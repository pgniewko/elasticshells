#ifndef VERTEXTRIANGLE_H
#define	VERTEXTRIANGLE_H

#include "Vertex.h"
#include "Triangle.h"

class Vertex; 

class VertexTriangle {
public:
    VertexTriangle();
    VertexTriangle(Vertex*, Vertex*, Vertex*);
    VertexTriangle(const VertexTriangle& orig);
    virtual ~VertexTriangle();
    void setId(int);
    double area();
    
    Vertex * a;
    Vertex * b;
    Vertex * c;
    int id;
};

#endif	/* VERTEXTRIANGLE_H */

