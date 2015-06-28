#include "TotalCellsArea.h"

TotalCellsArea::TotalCellsArea() {}

TotalCellsArea::TotalCellsArea(const TotalCellsArea& orig) {}

TotalCellsArea::~TotalCellsArea() {}

double TotalCellsArea::totalCellArea(std::vector<Cell>& cells)
{
    int cellsnumber = cells.size();
    double total_area = 0.0;

    for (int i = 0; i < cellsnumber; i++)
    {
        total_area += cells[i].calcSurfaceArea();
    }

    return total_area;
}