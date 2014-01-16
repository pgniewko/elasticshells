#ifndef CELL_H
#define	CELL_H

#include <iostream>
#include <cmath>
#include "Vector3D.h"

class Cell
{
    protected:
        double* params;

    public:
        Cell() : r(0, 0, 0), p(0, 0, 0), f(0, 0, 0), mass(0.0) {};
        Cell(double x, double y, double z);
//        Cell(const Cell& orig);
//        virtual ~Cell();

        virtual double calc_volume() {
            return 0.0;
        };

        Vector3D r;
        Vector3D p;
        Vector3D f;
        double mass;
};

#endif	/* CELL_H */

