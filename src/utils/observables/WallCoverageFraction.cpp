#include "WallCoverageFraction.h"

WallCoverageFraction::WallCoverageFraction(const char* name, const char* format) : Observer(name, format) {}

WallCoverageFraction::WallCoverageFraction(const WallCoverageFraction& orig) : Observer(orig) {}

WallCoverageFraction::~WallCoverageFraction() {}

void WallCoverageFraction::set_params(int num, ...)
{
    va_list arguments;
    va_start (arguments, num);
    d_param = va_arg(arguments, double);
    va_end( arguments );
}

void WallCoverageFraction::set_params(int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
}

double WallCoverageFraction::observe(Box& box, std::vector<Cell>& cells)
{
    double box_area = box.getArea(d_param);
    uint cells_number = cells.size();
    double coverage = 0.0;

    for (uint i = 0; i < cells_number; i++)
    {
        coverage += cells[i].contactArea(box, d_param);
    }

    coverage /= box_area;
    return coverage;
}

DerivedRegister<WallCoverageFraction> WallCoverageFraction::reg("WallCoverageFraction");