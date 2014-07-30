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
        Cell() : r(0, 0, 0), p(0, 0, 0), cf(0, 0, 0), mass(0.0) {};
        Cell(double x, double y, double z);
//        Cell(const Cell& orig);
//        virtual ~Cell();

        virtual double calc_volume() {
            return 0.0;
        };
        virtual double cforce_magnitude() {
            return cf.length();
        }
        virtual Vector3D cforce_direction() {
            return cf. direction();
        };
        virtual Vector3D reset_cforce();

        Vector3D r;
        Vector3D p;
        Vector3D cf;
        double mass;
};

#endif	/* CELL_H */
