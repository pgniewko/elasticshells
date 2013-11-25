/* 
 * File:   Cell.h
 * Author: pablo
 *
 * Created on November 25, 2013, 11:29 AM
 */

#ifndef CELL_H
#define	CELL_H

class Cell 
{

public:
    Cell();
    Cell(const Cell& orig);
    virtual ~Cell();
    void set_coor(double, double, double);
    void set_radius(double);
    double calc_volume();
private:
    double coor[3];
    double rot[3];
    double radius;
};

#endif	/* CELL_H */

