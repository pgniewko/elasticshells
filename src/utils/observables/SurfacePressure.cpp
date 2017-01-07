#include "SurfacePressure.h"

SurfacePressure::SurfacePressure(const char* name, const char* format) : Observer(name, format) {}

SurfacePressure::SurfacePressure(const SurfacePressure& orig) : Observer(orig) {}

SurfacePressure::~SurfacePressure() {}

void SurfacePressure::set_params(const int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
}

double SurfacePressure::observe(const Box& box, std::vector<Cell>& cells)
{
    if (box.pbc)
    {
        return 0.0;
    }
    
    double totalForce = SurfaceForce::calcTotalForce(box, cells);
    double area = box.getArea(d_param);
    return totalForce / area;
}

DerivedRegister<SurfacePressure> SurfacePressure::reg("SurfacePressure");