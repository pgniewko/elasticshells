#ifndef SIMPLETRIANGULATION_H
#define	SIMPLETRIANGULATION_H

#include "Triangulation.h"

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

