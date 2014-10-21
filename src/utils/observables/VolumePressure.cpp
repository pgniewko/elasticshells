#include "VolumePressure.h"

VolumePressure::VolumePressure() {}

VolumePressure::VolumePressure(const VolumePressure& orig) {}

VolumePressure::~VolumePressure() {}

double VolumePressure::calcPressure(Box& box, vector<Cell>& cells)
{
    double V = box.getVolume();
    double sab[] = {0, 0, 0, 0, 0, 0};
    int numofcells = cells.size();
    for (int i = 0; i < numofcells; i++)
    {
        cells[i].calcStressTensor(box, sab);
    }
    
    double pressure = -(sab[0] + sab[1] + sab[2]);
    pressure /= 3.0;
    pressure /= V;
    return pressure;
}