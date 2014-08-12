#ifndef VERTEX_H
#define	VERTEX_H

#define NEIGH_MAX 10
#define TRIAN_MAX 10

#include "Vector3D.h"
//#include "VertexTriangle.h"

class Vertex
{
public:
    Vertex();
    Vertex(double, double, double);
    Vertex(Vector3D);
    Vertex(const Vertex& orig);
    virtual ~Vertex();
    int setId(int);
    int getId();
    double setMass(double);
    double getMass();
    void printVertex();
    
    Vector3D force;
    Vector3D velocity;
    Vector3D xyz;
    
    //Vertex neighbors[NEIGH_MAX];
    //double k0 [NEIGH_MAX];
    //VertexTriangle vtris [TRIAN_MAX];
    
private:
    int id;
    double mass;
};

#endif	/* VERTEX_H */

