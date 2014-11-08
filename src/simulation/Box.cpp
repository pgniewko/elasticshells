#include "Box.h"

Box::Box(double bsx, double bsy, double bsz) : pbc(false), ecw(0.0),
    x(bsx), y(bsy), z(bsz), xs(bsx), ys(bsy), zs(bsz),
    xe(0.5 * bsx), ye(0.5 * bsy), ze(0.5 * bsz),
    dx(0), dy(0), dz(0)
{
}

Box::Box(double bsx, double bsy, double bsz, double dbs) : pbc(false), ecw(0.0),
    x(bsx), y(bsy), z(bsz), xs(bsx), ys(bsy), zs(bsz),
    xe(0.5 * bsx), ye(0.5 * bsy), ze(0.5 * bsz),
    dx(dbs), dy(dbs), dz(dbs)
{
}

Box::Box(const Box& orig) : pbc(orig.pbc), ecw(orig.ecw),
    x(orig.x), y(orig.y), z(orig.z), xs(orig.xs), ys(orig.ys), zs(orig.zs),
    xe(orig.xe), ye(orig.ye), ze(orig.ze),
    dx(orig.dx), dy(orig.dy), dz(orig.dz)
{
}

Box::~Box() {}

void Box::setPbc(bool pbcf)
{
    pbc = pbcf;
}

void Box::setEcw(double e)
{
    ecw = e;
}

void Box::resize()
{
    if (x + dx >= xe)
    {
        x += dx;
    }

    if (y + dy >= ye)
    {
        y += dy;
    }

    if (z + dz >= ze)
    {
        z += dz;
    }
}

double Box::getVolume()
{
    return 2.0 * x * 2.0 * y * 2.0 * z;
}

double Box::getArea()
{
    return 2 * (4 * x * y + 4 * x * z + 4 * y * z);
}

void Box::setX(const double newx)
{
    x = newx;
}

double Box::getX()
{
    return x;
}

void Box::setY(const double newy)
{
    y = newy;
}

double Box::getY()
{
    return y;
}

void Box::setZ(const double newz)
{
    z = newz;
}

double Box::getZ()
{
    return z;
}

void Box::setDx(const double newdx)
{
    dx = newdx;
}

double Box::getDx()
{
    return dx;
}

void Box::setDy(const double newdy)
{
    dy = newdy;
}

double Box::getDy()
{
    return dy;
}

void Box::setDz(const double newdz)
{
    dz = newdz;
}

double Box::getDz()
{
    return dz;
}

void Box::setXstart(const double xst)
{
    xs = xst;
}

void Box::setYstart(const double yst)
{
    ys = yst;
}

void Box::setZstart(const double zst)
{
    zs = zst;
}

void Box::setXend(const double xend)
{
    xe = xend;
}

void Box::setYend(const double yend)
{
    ye = yend;
}

void Box::setZend(const double zend)
{
    ze = zend;
}

double Box::getXstart()
{
    return xs;
}
double Box::getYstart()
{
    return ys;
}

double Box::getZstart()
{
    return zs;
}

double Box::getXend()
{
    return xe;
}
double Box::getYend()
{
    return ye;
}

double Box::getZend()
{
    return ze;
}