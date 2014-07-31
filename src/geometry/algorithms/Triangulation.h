#ifndef TRIANGULATION_H
#define	TRIANGULATION_H

#include <list>
#include "../Triangle.h"

using namespace std;

class Triangulation {
public:
    Triangulation();
    Triangulation(const Triangulation& orig);
    virtual ~Triangulation();
    virtual list<Triangle> triangulate() =0;
private:
    std::list<Triangle> tris;

};

#endif	/* TRIANGULATION_H */

