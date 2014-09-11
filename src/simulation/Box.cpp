#include "Box.h"

Box::Box(double bsx, double bsy, double bsz) : x(bsx), y(bsy), z(bsz), 
       dx(0), dy(0), dz(0), xe(0.5*bsx), ye(0.5*bsy), ze(0.5*bsz) {}

Box::Box(double bsx, double bsy, double bsz, double dbs) : x(bsx), y(bsy), z(bsz), 
        dx(dbs), dy(dbs), dz(dbs), xe(0.5*bsx), ye(0.5*bsy), ze(0.5*bsz) {}

Box::Box(const Box& orig) : x(orig.x), y(orig.y), z(orig.z), 
        xe(orig.xe), ye(orig.ye), ze(orig.ze), dx(orig.dx), dy(orig.dy), dz(orig.dz) {}

Box::~Box() {}

void Box::resize()
{
    //std::cout << "resized: " << " x= " << x << " y= "<< y << " z=" << z;
    //std::cout << " xe= " << xe << " ye= "<< ye << " ze=" << ze << std::endl;
    if (x + dx >= xe) {x += dx;}
    
    if (y + dy >= ye) {y += dy;}
    
    if (z + dz >= ze) {z += dz;}
    //std::cout << "resized: " << " x= " << x << " y= "<< y << " z=" << z << std::endl;
}

double Box::getVolume()
{
    return 2.0 * x * 2.0 * y * 2.0 * z;
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

void Box::setXend(const double xend) {xe = xend;}
void Box::setYend(const double yend) {ye = yend;}
void Box::setZend(const double zend) {ze = zend;}

double Box::getXend() {return xe;}
double Box::getYend() {return ye;}
double Box::getZend() {return ze;}

