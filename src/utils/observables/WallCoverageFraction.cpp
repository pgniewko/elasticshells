#include "WallCoverageFraction.h"

WallCoverageFraction::WallCoverageFraction(const char* name, const char* format) : Observer(name, format)
{}

WallCoverageFraction::WallCoverageFraction(const WallCoverageFraction& orig) : Observer(orig.observer_name, orig.output_format), rv(orig.rv)
{}

WallCoverageFraction::~WallCoverageFraction() 
{}

void SurfacePressure::set_params(int num, ...)
{
    va_list arguments;
    va_start (arguments, num);
    rv = va_arg(arguments, double);
    va_end( arguments );
}

void SurfacePressure::set_params(int num, std::vector<std::string> args_)
{
    rv = strtod(args_[ num + 0 ].c_str(), NULL);
}

double WallCoverageFraction::observe(Box& box, std::vector<Cell>& cells)
{
    double box_area = box.getArea(rv);
    int cells_number = cells.size();
    double coverage = 0.0;

    for (int i = 0; i < cells_number; i++)
    {
        coverage += cells[i].contactArea(box);
    }

    coverage /= box_area;
    return coverage;
}

DerivedRegister<WallCoverageFraction> WallCoverageFraction::reg("WallCoverageFraction");