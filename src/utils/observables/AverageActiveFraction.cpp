#include "AverageActiveFraction.h"

AverageActiveFraction::AverageActiveFraction(const char* name, const char* format) : Observer(name, format) {}

AverageActiveFraction::AverageActiveFraction(const AverageActiveFraction& orig) : Observer(orig) {}

AverageActiveFraction::~AverageActiveFraction() {}

void AverageActiveFraction::set_params(const int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
    return;
};

double AverageActiveFraction::observe(const Box& box, std::vector<Shell>& cells, const DomainList& dl)
{

    double total_active_f = 0.0;
    double active_a = 0.0;
    double cell_a = 0.0;
    double cell_f = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        cells[i].calcCM();
    }

    for (uint i = 0; i < cells.size(); i++)
    {
        active_a = cells[i].activeArea(box, cells, d_param);
        cell_a = cells[i].calcSurfaceArea(d_param);
        cell_f = active_a / cell_a;
        total_active_f += cell_f;
    }

    total_active_f /= cells.size();


    return total_active_f;
}

DerivedRegister<AverageActiveFraction> AverageActiveFraction::reg("AverageActiveFraction");
