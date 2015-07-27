#include "SurfacePressure.h"

SurfacePressure::SurfacePressure(const char* name, const char* format) : Observer(name, format)
{
//    std::cout << "SurfacePressure:Buduja mnie"<<std::endl;
//    DerivedRegister<SurfacePressure> SurfacePressure::reg("SurfacePressure");
//    rv=0.10;
//    std::cout << "rv=" << rv << std::endl;
}

SurfacePressure::SurfacePressure(const SurfacePressure& orig) : Observer(orig.observer_name, orig.output_format)
{std::cout << "SurfacePressure:Kopjuja mnie"<<std::endl;}

SurfacePressure::~SurfacePressure() 
{std::cout << "SurfacePressure:Niszcza mnie"<<std::endl;}

double SurfacePressure::observe(Box& box, std::vector<Cell>& cells)
{
    double totalForce = SurfaceForce::calcForces(box, cells) ;
    double area = box.getArea(rv);
    //std::cout << "Moj rv=" << rv << std::endl;
    return totalForce / area;
}

void SurfacePressure::set_params(int num, ...)
{
    va_list arguments;
    va_start (arguments, num);
    rv = va_arg(arguments, double);
    va_end( arguments );
//    std::cout << "Moj rv=" << rv << std::endl;
}

DerivedRegister<SurfacePressure> SurfacePressure::reg("SurfacePressure");