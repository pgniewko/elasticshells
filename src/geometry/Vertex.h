#ifndef VERTEX_H
#define	VERTEX_H

#define NEIGH_SIZE 10

#include "Vector3D.h"

class Vertex
{
public:
//    Vertex();
    Vertex(double, double, double);
    Vertex(Vector3D);
    Vertex(const Vertex& orig);
    virtual ~Vertex();
    int setId(int);
    int getId();
    double setMass(double);
    double getMass();
    void addNeigh(int);
    void addTriangle(int);
    void printVertex();
    
    Vector3D force;
    Vector3D velocity;
    Vector3D xyz;
    int neigsIds [NEIGH_SIZE];
    double k0 [NEIGH_SIZE];
    int trisIds [NEIGH_SIZE];
    int bondedNo;
    int noTris;
//    int nbneigs [NEIGH_SIZE];
//    int nbNo;
//    int bNo;
//    double mass;
private:
    int id;
    double mass;
};

#endif	/* VERTEX_H */

