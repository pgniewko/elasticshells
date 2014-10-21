#include "VolumeFraction.h"

VolumeFraction::VolumeFraction() {}

VolumeFraction::VolumeFraction(const VolumeFraction& orig) {}

VolumeFraction::~VolumeFraction() {}

double VolumeFraction::caclVolumeFraction(Box& box, vector<Cell>& cells)
{
    double boxVolume = box.getVolume();
    double cellsVolume = 0.0;
    int numofcells = cells.size();
    for (int i = 0; i < numofcells; i++)
    {
        cellsVolume += cells[i].calcVolume();
    }
    double volumeFraction = cellsVolume / boxVolume;
    return volumeFraction;
}

double VolumeFraction::caclCellsVolume(vector<Cell>& cells)
{

    double cellsVolume = 0.0;
    int numofcells = cells.size();
    for (int i = 0; i < numofcells; i++)
    {
        cellsVolume += cells[i].calcVolume();
    }
    return cellsVolume;
}
