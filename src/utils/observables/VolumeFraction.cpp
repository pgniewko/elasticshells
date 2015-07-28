#include "VolumeFraction.h"

VolumeFraction::VolumeFraction(const char* name, const char* format) : Observer(name, format)
{}

VolumeFraction::VolumeFraction(const VolumeFraction& orig) : Observer(orig.observer_name, orig.output_format), rv(orig.rv)
{}

VolumeFraction::~VolumeFraction() {}

double VolumeFraction::observe(Box& box, std::vector<Cell>& cells)
{
    double boxVolume = box.getVolume(0.0);
    double cellsVolume = calcCellsVolume(cells);
    return (cellsVolume / boxVolume);
}

double VolumeFraction::calcCellsVolume(std::vector<Cell>& cells)
{
    double cellsVolume = 0.0;
    int numofcells = cells.size();

    for (int i = 0; i < numofcells; i++)
    {
        cellsVolume += cells[i].calcVolume(rv);
    }

    return cellsVolume;
}

void VolumeFraction::set_params(int num, ...)
{
    va_list arguments;
    va_start (arguments, num);
    rv = va_arg(arguments, double);
    va_end( arguments );
}

void VolumeFraction::set_params(int num, std::vector<std::string> args_)
{
    rv = strtod(args_[ num + 0 ].c_str(), NULL);
}

DerivedRegister<VolumeFraction> VolumeFraction::reg("VolumeFraction");
