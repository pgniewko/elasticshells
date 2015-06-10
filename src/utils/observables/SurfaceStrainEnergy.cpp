#include "SurfaceStrainEnergy.h"

SurfaceStrainEnergy::SurfaceStrainEnergy() {}

SurfaceStrainEnergy::SurfaceStrainEnergy(const SurfaceStrainEnergy& orig) {}

SurfaceStrainEnergy::~SurfaceStrainEnergy() {}

double SurfaceStrainEnergy::calcSurfaceEnergy(std::vector<Cell>& cells)
{
    double surf_energy = 0.0;
    for (int i = 0; i < cells.size();)
    {
        surf_energy += cells[i].surfaceStrainEnergy();
    }
    return surf_energy;
}