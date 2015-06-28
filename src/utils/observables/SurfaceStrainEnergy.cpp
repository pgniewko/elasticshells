#include "SurfaceStrainEnergy.h"

SurfaceStrainEnergy::SurfaceStrainEnergy() {}

SurfaceStrainEnergy::SurfaceStrainEnergy(const SurfaceStrainEnergy& orig) {}

SurfaceStrainEnergy::~SurfaceStrainEnergy() {}

double SurfaceStrainEnergy::calcSurfaceEnergy(std::vector<Cell>& cells)
{
    double surf_energy = 0.0;
    int cellsnumber = cells.size();

    for (int i = 0; i < cellsnumber; i++)
    {
        surf_energy += cells[i].surfaceStrainEnergy();
    }

    return surf_energy;
}