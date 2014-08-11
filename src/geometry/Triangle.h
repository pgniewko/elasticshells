#ifndef TRIANGLE_H
#define	TRIANGLE_H

#include "Vector3D.h"

class Triangle {
public:
//    Triangle();
//    Triangle(Vector3D, Vector3D, Vector3D);
    Triangle(Vector3D&, Vector3D&, Vector3D&);
    Triangle(const Triangle& orig);
    virtual ~Triangle();
    double area() const;
    void setId(int);
    int getId();
    void printTriangle();
    Vector3D * a, * b, * c;
    int id;
     
private:
   

};

#endif	/* TRIANGLE_H */

