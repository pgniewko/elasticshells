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

void Cell::set_rot(double a1, double a2, double a3){
    rot[0] = a1;
    rot[1] = a2;
    rot[2] = a3;
}

Cell::Cell() {
}

Cell::Cell(const Cell& orig) {
}

Cell::~Cell() {
}

