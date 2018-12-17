#include "AverageActiveArea.h"

AverageActiveArea::AverageActiveArea(const char* name, const char* format) : Observer(name, format) {}

AverageActiveArea::AverageActiveArea(const AverageActiveArea& orig) : Observer(orig) {}

AverageActiveArea::~AverageActiveArea() {}

void AverageActiveArea::set_params(const int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
    return;
};

double AverageActiveArea::observe(const Box& box, std::vector<Shell>& shells, const DomainList& dl)
{
    double total_active_a = 0.0;

    for (uint i = 0; i < shells.size(); i++)
    {
        shells[i].calcCM();
    }

    for (uint i = 0; i < shells.size(); i++)
    {
        total_active_a += shells[i].activeArea(box, shells, d_param);
    }

    total_active_a /= shells.size();

    return total_active_a;
}

DerivedRegister<AverageActiveArea> AverageActiveArea::reg("AverageActiveArea");