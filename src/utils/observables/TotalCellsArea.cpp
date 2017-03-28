#include "TotalCellsArea.h"

TotalCellsArea::TotalCellsArea(const char* name, const char* format) : Observer(name, format) {}

TotalCellsArea::TotalCellsArea(const TotalCellsArea& orig) : Observer(orig) {}

TotalCellsArea::~TotalCellsArea() {}

void TotalCellsArea::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double TotalCellsArea::observe(const Box& box, std::vector<Cell>& cells, const DomainList& dl)
{
    uint cellsnumber = cells.size();
    double total_area = 0.0;

    for (uint i = 0; i < cellsnumber; i++)
    {
        total_area += cells[i].calcSurfaceArea();
    }

    return total_area;
}

DerivedRegister<TotalCellsArea> TotalCellsArea::reg("TotalCellsArea");