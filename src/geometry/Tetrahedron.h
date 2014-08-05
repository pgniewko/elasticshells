#include "Vector3D.h"

#ifndef TETRAHEDRON_H
#define	TETRAHEDRON_H

class Tetrahedron {
public:
//    Tetrahedron();
    Tetrahedron(Vector3D, Vector3D, Vector3D, Vector3D);
    Tetrahedron(const Tetrahedron& orig);
    virtual ~Tetrahedron();
    double volume() const;
    Vector3D a, b, c, d;
private:

};

#endif	/* TETRAHEDRON_H */

