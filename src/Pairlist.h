#ifndef PAIRLIST_H
#define PAIRLIST_H

#include <vector>
#include <iostream>
#include "Vector3D.h"
#include "Distance.h"

using namespace std;

typedef vector<int> IntList;
typedef vector<IntList> BoxList;

class Pairlist {
public:
    Pairlist(double mindist, int steps = 1);
    virtual ~Pairlist();
    
    void compute( const  Distance& dist, const Vector3D* pos, int  np, int* aindex = 0);
    int findbox(Vector3D p, const Distance& dist ) const;
    int num_boxes() const; // const {return xb * yb * zb;}
    BoxList boxlist;
    BoxList nbrlist;
    BoxList fullnbrlist;
    
private:
    int n;
    double mindist;
    Vector3D margin;
    Vector3D lmin;
    Vector3D pairdist;
    int steps;
    int xb, yb, zb, xytotb;
    void set_neighbors(const Distance& dist);
};
#endif