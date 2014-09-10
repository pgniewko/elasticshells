#ifndef SIMPLETRIANGULATION_H
#define	SIMPLETRIANGULATION_H

#include <iostream>    /* cout, cin */

#include "Environment.h"
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
    list<Triangle> triangulate(double);
private:
    void createCube();
    void subdivide();
    int depth;

};

#endif	/* SIMPLETRIANGULATION_H */
