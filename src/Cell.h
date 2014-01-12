#ifndef CELL_H
#define	CELL_H

#include <iostream>
#include <cmath>

class Cell {
    protected:
        double coor[3];
        double rot[3];
        double* params;

    public:
        Cell();
        Cell(const Cell& orig);
        virtual ~Cell();

        virtual void set_coor(double, double, double);
        virtual void set_rot(double, double, double);
        virtual void set_params(double*, int) = 0;
        virtual double calc_volume() = 0;
};

#endif	/* CELL_H */

