#include "VolumeFraction.h"

VolumeFraction::VolumeFraction() {}

VolumeFraction::VolumeFraction(const VolumeFraction& orig) {}

VolumeFraction::~VolumeFraction() {}

double VolumeFraction::caclVolumeFraction(Box& box, std::vector<Cell>& cells, double rv)
{
    //double boxVolume = box.getVolume(rv);
    double boxVolume = box.getVolume(0.0);
    double cellsVolume = caclCellsVolume(cells, rv);
    return (cellsVolume / boxVolume);
}

double VolumeFraction::caclCellsVolume(std::vector<Cell>& cells, double eps)
{
    double cellsVolume = 0.0;
    int numofcells = cells.size();

    for (int i = 0; i < numofcells; i++)
    {
        cellsVolume += cells[i].calcVolume(eps);
    }

    return cellsVolume;
}
