#include "TotalContactArea.h"

TotalContactArea::TotalContactArea(const char* name, const char* format) : Observer(name, format) {}

TotalContactArea::TotalContactArea(const TotalContactArea& orig) : Observer(orig) {}

TotalContactArea::~TotalContactArea() {}

void TotalContactArea::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double TotalContactArea::observe(const Box& box, std::vector<Shell>& shells, const DomainList& dl)
{
    uint shellsnumber = shells.size();
    double total_contact_a = 0.0;

    for (uint i = 0; i < shellsnumber - 1; i++)
    {
        for (uint j = i + 1; j < shellsnumber; j++)
        {
            total_contact_a += shells[i].contactArea(shells[j], box);
        }
    }

    return total_contact_a;
}

DerivedRegister<TotalContactArea> TotalContactArea::reg("TotalContactArea");