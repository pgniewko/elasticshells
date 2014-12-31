#include "VolumeFraction.h"

VolumeFraction::VolumeFraction() {}

VolumeFraction::VolumeFraction(const VolumeFraction& orig) {}

VolumeFraction::~VolumeFraction() {}

double VolumeFraction::caclVolumeFraction(Box& box, std::vector<Cell>& cells, double rv)
{
    double boxVolume = box.getVolume(rv);
    double cellsVolume = caclCellsVolume(cells);
    return (cellsVolume / boxVolume);
}

double VolumeFraction::caclCellsVolume(std::vector<Cell>& cells)
{
    double cellsVolume = 0.0;
    int numofcells = cells.size();

    for (int i = 0; i < numofcells; i++)
    {
        cellsVolume += cells[i].calcVolume();
    }

    return cellsVolume;
}
