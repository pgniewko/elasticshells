/* 
 * File:   Cell.cpp
 * Author: pablo
 * 
 * Created on November 25, 2013, 11:29 AM
 */

#include "Cell.h"

void Cell::set_coor(double x, double y, double z)
{
    coor[0] = x;
    coor[1] = y;
    coor[2] = z;
}

void Cell::set_radius(double r)
{
    radius = r;
}

double Cell::calc_volume(){
    return (4.0*3.14 * radius*radius*radius/3);
}


Cell::Cell() {
}

Cell::Cell(const Cell& orig) {
}

Cell::~Cell() {
}

