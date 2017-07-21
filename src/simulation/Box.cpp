#include "Box.h"

utils::Logger Box::box_logger("box_logger");

Box::Box(double bsx, double bsy, double bsz) : pbc(false),
    x(bsx), y(bsy), z(bsz), x_min(bsx), y_min(bsy), z_min(bsz),
    E_box(0.0), nu(0.0) {}

Box::Box(const Box& orig) : pbc(orig.pbc),
    x(orig.x), y(orig.y), z(orig.z),
    x_min(orig.x_min), y_min(orig.y_min), z_min(orig.z_min),
    E_box(orig.E_box), nu(orig.nu), my_schedule(orig.my_schedule) {}

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

bool Box::resize(double vf_)
{
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;

    my_schedule.execute(dx, dy, dz, vf_);

    if (x + dx >= x_min)
    {
        x += dx;
    }

    if (y + dy >= y_min)
    {
        y += dy;
    }

    if (z + dz >= z_min)
    {
        z += dz;
    }

    if (dx == 0 && dy == 0 && dz == 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

double Box::getVolume(const double rv) const
{
    return 2.0 * (x - rv) * 2.0 * (y - rv) * 2.0 * (z - rv);
}

double Box::getArea(const double rv) const
{
    double xc = x - rv;
    double yc = y - rv;
    double zc = z - rv;
    return 2 * (4 * xc * yc + 4 * xc * zc + 4 * yc * zc);
}

void Box::setX(const double newx)
{
    x = newx;
}

double Box::getX() const
{
    return x;
}

void Box::setY(const double newy)
{
    y = newy;
}

double Box::getY() const
{
    return y;
}

void Box::setZ(const double newz)
{
    z = newz;
}

double Box::getZ() const
{
    return z;
}

void Box::setXmin(const double xend)
{
    x_min = xend;
}

void Box::setYmin(const double yend)
{
    y_min = yend;
}

void Box::setZmin(const double zend)
{
    z_min = zend;
}

double Box::getXmin() const
{
    return x_min;
}
double Box::getYmin() const
{
    return y_min;
}

double Box::getZmin() const
{
    return z_min;
}

double Box::getXEdge(const double rv) const
{
    return 2 * (x - rv);
}

double Box::getYEdge(const double rv) const
{
    return 2 * (y - rv);
}

double Box::getZEdge(const double rv) const
{
    return 2 * (z - rv);
}

double Box::getE() const
{
    return E_box;
}

double Box::getNu() const
{
    return nu;
}

void Box::configureScheduler(std::string schf)
{
    my_schedule.setFileName(schf);
    my_schedule.registerSchedules();
    my_schedule.configureSchedule();
}

void Box::setDefaultSchedule(int ns, int in, double _dx, double _dy, double _dz, double _rx, double _ry, double _rz)
{
    my_schedule.setDefault(ns, in, _dx, _dy, _dz, _rx, _ry, _rz);
}

void Box::saveRemainingSchedule()
{
    my_schedule.saveRemainingSchedule();
}

bool Box::nthTodo()
{
    return my_schedule.nthTodo();
}