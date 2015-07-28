#include "SurfacePressure.h"

SurfacePressure::SurfacePressure(const char* name, const char* format) : Observer(name, format)
{}

SurfacePressure::SurfacePressure(const SurfacePressure& orig) : Observer(orig)
{}

SurfacePressure::~SurfacePressure()
{}

void SurfacePressure::set_params(int num, ...)
{
    va_list arguments;
    va_start (arguments, num);
    d_param = va_arg(arguments, double);
    va_end( arguments );
}

void SurfacePressure::set_params(int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
}

double SurfacePressure::observe(Box& box, std::vector<Cell>& cells)
{
    double totalForce = SurfaceForce::calcTotalForce(box, cells);
    double area = box.getArea(d_param);
    return totalForce / area;
}

DerivedRegister<SurfacePressure> SurfacePressure::reg("SurfacePressure");