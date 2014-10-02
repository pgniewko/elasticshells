#include "SurfacePressure.h"

SurfacePressure::SurfacePressure() {}

SurfacePressure::SurfacePressure(const SurfacePressure& orig) {}

SurfacePressure::~SurfacePressure() {}

double SurfacePressure::calcPressure(Box& box, vector<Cell>& cells)
{
    double area = box.getArea();
    double forcemagn = 0.0;
    int numofcells = cells.size();
    for (int i = 0; i < numofcells; i++)
    {
        forcemagn += cells[i].calcBoxForces(box);
    }
    double pressure = forcemagn / area;
    return pressure;
}
