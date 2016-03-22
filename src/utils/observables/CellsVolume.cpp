#include "CellsVolume.h"

CellsVolume::CellsVolume (const char* name, const char* format) : Observer(name, format) {}

CellsVolume::CellsVolume(const CellsVolume& orig) : Observer(orig) {}

CellsVolume::~CellsVolume() {}

void CellsVolume::set_params(const int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
}

double CellsVolume::observe(const Box& box, std::vector<Cell>& cells)
{
    return calcCellsVolume(cells);
}

double CellsVolume::calcCellsVolume(std::vector<Cell>& cells)
{
    double cellsVolume = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        cellsVolume += cells[i].calcVolume(d_param);
    }

    return cellsVolume;
}

DerivedRegister<CellsVolume> CellsVolume::reg("CellsVolume");