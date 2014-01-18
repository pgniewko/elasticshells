#ifndef FORCE_H
#define	FORCE_H

#include "Vector3D.h"
#include "Cell.h"
#include <iostream>
using namespace std;

class Force {
public:
    Force();
    Force(const Force& orig);
    virtual ~Force();
    void print() {cout << "print func" <<endl;};
    virtual Vector3D eval(Vector3D r_kl, Cell cell_k, Cell cell_l) =0;
    inline virtual double magn(Vector3D r_kl, Cell cell_k, Cell cell_l) =0;
private:

};

#endif	/* FORCE_H */

