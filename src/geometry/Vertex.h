#ifndef VERTEX_H
#define	VERTEX_H

#define NEIGH_MAX 10
#define TRIAN_MAX 10

#include "Vector3D.h"

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
    
    void addNeighbor(int, double);
    void addTriangle(int);
    void voidForce();
    void voidVelocity();
    
    bool isBonded(int);
    
    Vector3D force;
    Vector3D velocity;
    Vector3D xyz;
    
    int neighbors[NEIGH_MAX];
    int vertextri[TRIAN_MAX];
    double R0[NEIGH_MAX];
    
    int nneigh;
    int ntrian;

private:
    int id;
    double mass;
};

#endif	/* VERTEX_H */

