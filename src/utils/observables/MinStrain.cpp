#include "MinStrain.h"

MinStrain::MinStrain (const char* name, const char* format) : Observer(name, format) {}

MinStrain::MinStrain (const MinStrain& orig) : Observer(orig) {}

MinStrain::~MinStrain () {}

void MinStrain::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double MinStrain::observe(const Box& box, std::vector<Cell>& cells)
{
    double min_strain = 10.0;
    uint cellsnumber = cells.size();

    for (uint i = 0; i < cellsnumber; i++)
    {
        min_strain = std::min(min_strain, cells[i].minStrain());
    }

    return min_strain;
}

DerivedRegister<MinStrain> MinStrain::reg("MinStrain");
