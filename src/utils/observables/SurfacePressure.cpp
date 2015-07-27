#include "SurfacePressure.h"

SurfacePressure::SurfacePressure(const char* name, const char* format) : Observer(name, format)
{
}

SurfacePressure::SurfacePressure(const SurfacePressure& orig) : Observer(orig.observer_name, orig.output_format)
{
}

SurfacePressure::~SurfacePressure() 
{
}

double SurfacePressure::observe(Box& box, std::vector<Cell>& cells)
{
    double totalForce = SurfaceForce::calcForces(box, cells) ;
    double area = box.getArea(rv);
    return totalForce / area;
}

void SurfacePressure::set_params(int num, ...)
{
    va_list arguments;
    va_start (arguments, num);
    rv = va_arg(arguments, double);
    va_end( arguments );
}

DerivedRegister<SurfacePressure> SurfacePressure::reg("SurfacePressure");