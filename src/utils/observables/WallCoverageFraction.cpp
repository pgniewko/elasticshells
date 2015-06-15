#include "WallCoverageFraction.h"

WallCoverageFraction::WallCoverageFraction() {}

WallCoverageFraction::WallCoverageFraction(const WallCoverageFraction& orig) {}

WallCoverageFraction::~WallCoverageFraction() {}

double WallCoverageFraction::wallsCoverage(Box& box, std::vector<Cell>& cells, double rv)
{
    double box_area = box.getArea(rv);
    //double box_area = box.getArea(0.0);
    int cells_number = cells.size();
    double coverage = 0.0;
    for (int i = 0; i < cells_number; i++)
    {
        coverage += cells[i].contactArea(box);
    }
    coverage /= box_area;
    return coverage;
}