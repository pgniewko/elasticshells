#ifndef TRIANGULATION_H
#define	TRIANGULATION_H

#include <list>
#include <iostream> 
#include "../Triangle.h"

using namespace std;

class Triangulation {
public:
    Triangulation();
    Triangulation(const Triangulation& orig);
    virtual ~Triangulation();
    virtual list<Triangle> triangulate() =0;
    void saveTriangulatedSurface(const char*, bool);
    void saveRenderingScript(const char*, const char*);
protected:
    list<Triangle> tris;

};

#endif	/* TRIANGULATION_H */
