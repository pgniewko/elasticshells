#ifndef VERTEX_H
#define	VERTEX_H

#include <stdlib.h>

#include "Environment.h"
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
    bool isNeighbor(int);
    void addTriangle(int);
    
    void voidForce();
    void voidVelocity();
    
    Vector3D xyz;
    Vector3D force;
    Vector3D velocity;
    
    Vector3D tmp_xyz;
    Vector3D tmp_force;
    Vector3D tmp_velocity;
    
    int neighbors[NEIGH_MAX];
    int vertextri[TRIAN_MAX];
    double R0[NEIGH_MAX];
    
    int nbvertices[NBNEI_MAX];
    int nbcellid[NBNEI_MAX];
    
    int nneigh;
    int ntrian;
    int nbneigh;

private:
    int id;
    double mass;
};

#endif	/* VERTEX_H */
