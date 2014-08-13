#ifndef SIMPLETRIANGULATION_H
#define	SIMPLETRIANGULATION_H

#include <iostream>    /* cout, cin */

#include "Triangulation.h"

using namespace std;

class SimpleTriangulation : public Triangulation
{
public:
    SimpleTriangulation();
    SimpleTriangulation(int);
    SimpleTriangulation(const SimpleTriangulation& orig);
    virtual ~SimpleTriangulation();
    list<Triangle> triangulate();
private:
    void createCube();
    void subdivide();
    int depth;

};

#endif	/* SIMPLETRIANGULATION_H */
