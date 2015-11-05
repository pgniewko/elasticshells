#include "TotalStrainEnergy.h"

TotalStrainEnergy::TotalStrainEnergy(const char* name, const char* format) : Observer(name, format) {}

TotalStrainEnergy::TotalStrainEnergy(const TotalStrainEnergy& orig) : Observer(orig) {}

TotalStrainEnergy::~TotalStrainEnergy() {}

void TotalStrainEnergy::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double TotalStrainEnergy::observe(const Box& box, std::vector<Cell>& cells)
{
    double strain_energy = 0.0;
    uint cellsnumber = cells.size();

    for (uint i = 0; i < cellsnumber; i++)
    {
        strain_energy += cells[i].strainEnergy(box);
    }

    return strain_energy;
}

DerivedRegister<TotalStrainEnergy> TotalStrainEnergy::reg("TotalStrainEnergy");