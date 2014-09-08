#include "Box.h"

Box::Box(double bsx, double bsy, double bsz) : x(bsx), y(bsy), z(bsz), 
        dx(0), dy(0), dz(0) {}

Box::Box(double bsx, double bsy, double bsz, double dbs) : x(bsx), y(bsy), z(bsz), 
        dx(dbs), dy(dbs), dz(dbs) {}

Box::Box(const Box& orig) : x(orig.x), y(orig.y), z(orig.z), 
        dx(orig.dx), dy(orig.dy), dz(orig.dz){}

Box::~Box() {}

void Box::resize()
{
    x += dx;
    y += dy;
    z += dz;
}

void Box::setX(const double newx) {x = newx;}
double Box::getX() {return x;}
void Box::setY(const double newy) {y = newy;}
double Box::getY() {return y;}
void Box::setZ(const double newz) {z = newz;}
double Box::getZ() {return z;}

void Box::setDx(const double newdx) {dx = newdx;}
double Box::getDx() {return dx;}
void Box::setDy(const double newdy) {dy = newdy;}
double Box::getDy() {return dy;}
void Box::setDz(const double newdz) {dz = newdz;}
double Box::getDz() {return dz;}
