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

void Box::set_pbc(bool pbcf)
{
    pbc = pbcf;
}

void Box::set_E(double e)
{
    E_box = e;
}

void Box::set_nu(double n)
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

double Box::get_volume(const double rv) const
{
    return 2.0 * (x - rv) * 2.0 * (y - rv) * 2.0 * (z - rv);
}

double Box::get_area(const double rv) const
{
    double xc = x - rv;
    double yc = y - rv;
    double zc = z - rv;
    return 2 * (4 * xc * yc + 4 * xc * zc + 4 * yc * zc);
}

void Box::set_x(const double newx)
{
    x = newx;
}

double Box::get_x() const
{
    return x;
}

void Box::set_y(const double newy)
{
    y = newy;
}

double Box::get_y() const
{
    return y;
}

void Box::set_z(const double newz)
{
    z = newz;
}

double Box::get_z() const
{
    return z;
}

void Box::set_x_min(const double xend)
{
    x_min = xend;
}

void Box::set_y_min(const double yend)
{
    y_min = yend;
}

void Box::set_z_min(const double zend)
{
    z_min = zend;
}

double Box::get_x_min() const
{
    return x_min;
}
double Box::get_y_min() const
{
    return y_min;
}

double Box::get_z_min() const
{
    return z_min;
}

//double Box::getXEdge(const double rv) const
//{
//    return 2 * (x - rv);
//}
//
//double Box::getYEdge(const double rv) const
//{
//    return 2 * (y - rv);
//}
//
//double Box::getZEdge(const double rv) const
//{
//    return 2 * (z - rv);
//}

double Box::get_E() const
{
    return E_box;
}

double Box::get_nu() const
{
    return nu;
}

void Box::configure_scheduler(std::string schf)
{
    my_schedule.setFileName(schf);
    my_schedule.registerSchedules();
    my_schedule.configureSchedule();
}

void Box::set_default_schedule(int ns, int in, double _dx, double _dy, double _dz, double _rx, double _ry, double _rz)
{
    my_schedule.setDefault(ns, in, _dx, _dy, _dz, _rx, _ry, _rz);
}

void Box::save_remaining_schedule()
{
    my_schedule.saveRemainingSchedule();
}

bool Box::nthTodo()
{
    return my_schedule.nthTodo();
}