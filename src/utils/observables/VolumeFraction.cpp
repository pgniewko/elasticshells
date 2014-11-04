#include "VolumeFraction.h"

VolumeFraction::VolumeFraction() {}

VolumeFraction::VolumeFraction(const VolumeFraction& orig) {}

VolumeFraction::~VolumeFraction() {}

double VolumeFraction::caclVolumeFraction(Box& box, std::vector<Cell>& cells)
{
    double maxRbc = 0.0;
    
    for (int i = 0 ; i < cells.size(); i++)
    {
        maxRbc = std::max(maxRbc, cells[i].getRbc());
    }
    
    //double boxVolume = box.getVolume(maxRbc);
    double boxVolume = box.getVolume();
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
