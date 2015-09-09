#include "SurfaceStrainEnergy.h"

SurfaceStrainEnergy::SurfaceStrainEnergy(const char* name, const char* format) : Observer(name, format) {}

SurfaceStrainEnergy::SurfaceStrainEnergy(const SurfaceStrainEnergy& orig) : Observer(orig) {}

SurfaceStrainEnergy::~SurfaceStrainEnergy() {}

void SurfaceStrainEnergy::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double SurfaceStrainEnergy::observe(const Box& box, std::vector<Cell>& cells)
{
    double surf_energy = 0.0;
    int cellsnumber = cells.size();

    for (int i = 0; i < cellsnumber; i++)
    {
        surf_energy += cells[i].surfaceStrainEnergy();
    }

    return surf_energy;
}

DerivedRegister<SurfaceStrainEnergy> SurfaceStrainEnergy::reg("SurfaceStrainEnergy");