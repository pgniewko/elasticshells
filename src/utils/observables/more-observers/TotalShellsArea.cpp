#include "TotalShellsArea.h"

TotalShellsArea::TotalShellsArea(const char* name, const char* format) : Observer(name, format) {}

TotalShellsArea::TotalShellsArea(const TotalShellsArea& orig) : Observer(orig) {}

TotalShellsArea::~TotalShellsArea() {}

void TotalShellsArea::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double TotalShellsArea::observe(const Box& box, std::vector<Shell>& shells, const DomainList& dl)
{
    uint shellsrumber = shells.size();
    double total_area = 0.0;

    for (uint i = 0; i < shellsrumber; i++)
    {
        total_area += shells[i].calcSurfaceArea();
    }

    return total_area;
}

DerivedRegister<TotalShellsArea> TotalShellsArea::reg("TotalShellsArea");