#include "MaxStrain.h"

MaxStrain::MaxStrain (const char* name, const char* format) : Observer(name, format) {}

MaxStrain::MaxStrain (const MaxStrain& orig) : Observer(orig) {}

MaxStrain::~MaxStrain () {}

void MaxStrain::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double MaxStrain::observe(const Box& box, std::vector<Cell>& cells)
{
    double max_strain = -10.0;
    uint cellsnumber = cells.size();

    for (uint i = 0; i < cellsnumber; i++)
    {
        max_strain = std::max(max_strain, cells[i].maxStrain());
    }

    return max_strain;
}

DerivedRegister<MaxStrain> MaxStrain::reg("MaxStrain");