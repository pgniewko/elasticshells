#ifndef VERTEX_H
#define	VERTEX_H

#include "Vector3D.h"

class Vertex {
public:
    Vertex();
    Vertex(const Vertex& orig);
    virtual ~Vertex();
private:
    Vector3D force;
    Vector3D velocity;
    Vector3D xyz;
    int id;
    int bneigs [10];
    int nbneigs [10];
    double mass;

};

#endif	/* VERTEX_H */

