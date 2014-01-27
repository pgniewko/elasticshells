#ifndef FORCE_H
#define	FORCE_H

#include <iostream>

#include "Cell.h"
#include "Vector3D.h"

using namespace std;

class Force
{
    public:
        Force();
        Force(const Force& orig);
        virtual ~Force();
        virtual Vector3D eval(Vector3D r_kl, Cell cell_k, Cell cell_l) = 0;
        virtual double magn(Vector3D r_kl, Cell cell_k, Cell cell_l) = 0;
    private:

};

#endif	/* FORCE_H */

