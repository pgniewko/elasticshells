#include "VolumeFraction.h"

VolumeFraction::VolumeFraction(const char* name, const char* format) : Observer(name, format)
{}

VolumeFraction::VolumeFraction(const VolumeFraction& orig) : Observer(orig)
{}

VolumeFraction::~VolumeFraction() 
{}

void VolumeFraction::set_params(int num, ...)
{
    va_list arguments;
    va_start (arguments, num);
    d_param = va_arg(arguments, double);
    va_end( arguments );
}

void VolumeFraction::set_params(int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
}

double VolumeFraction::observe(Box& box, std::vector<Cell>& cells)
{
    double boxVolume = box.getVolume(0.0);
    double cellsVolume = calcCellsVolume(cells);
    return (cellsVolume / boxVolume);
}

double VolumeFraction::calcCellsVolume(std::vector<Cell>& cells)
{
    double cellsVolume = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        cellsVolume += cells[i].calcVolume(d_param);
    }

    return cellsVolume;
}

DerivedRegister<VolumeFraction> VolumeFraction::reg("VolumeFraction");