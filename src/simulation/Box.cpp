#include "Box.h"

Box::Box(double bsx, double bsy, double bsz) : pbc(false),
    x(bsx), y(bsy), z(bsz), xs(bsx), ys(bsy), zs(bsz),
    xe(bsx), ye(bsy), ze(bsz), dx(0), dy(0), dz(0), E_box(0.0), nu(0.0)
{
}

Box::Box(double bsx, double bsy, double bsz, double dbs) : pbc(false),
    x(bsx), y(bsy), z(bsz), xs(bsx), ys(bsy), zs(bsz),
    xe(bsx), ye(bsy), ze(bsz), dx(dbs), dy(dbs), dz(dbs), E_box(0.0), nu(0.0)
{
}

Box::Box(const Box& orig) : pbc(orig.pbc),
    x(orig.x), y(orig.y), z(orig.z), xs(orig.xs), ys(orig.ys), zs(orig.zs),
    xe(orig.xe), ye(orig.ye), ze(orig.ze),
    dx(orig.dx), dy(orig.dy), dz(orig.dz), E_box(orig.E_box), nu(orig.nu)
{
}

Box::~Box() {}

void Box::setPbc(bool pbcf)
{
    pbc = pbcf;
}

void Box::setEwall(double e)
{
    E_box = e;
}

void Box::setNu(double n)
{
    nu = n;
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

double Box::getVolume(double rv)
{
    return 2.0 * (x - rv) * 2.0 * (y - rv) * 2.0 * (z - rv);
}

double Box::getArea(double rv)
{
    double xc = x - rv;
    double yc = y - rv;
    double zc = z - rv;
    return 2 * (4 * xc * yc + 4 * xc * zc + 4 * yc * zc);
}

//double Box::getVolume()
//{
//    return 2.0 * x * 2.0 * y * 2.0 * z;
//}

//double Box::getArea()
//{
//    return 2 * (4 * x * y + 4 * x * z + 4 * y * z);
//}

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

double Box::getXEdge(double rv)
{
    return 2 * (x - rv);
}

double Box::getYEdge(double rv)
{
    return 2 * (y - rv);
}

double Box::getZEdge(double rv)
{
    return 2 * (z - rv);
}

double Box::getE()
{
    return E_box;
}

double Box::getNu()
{
    return nu;
}