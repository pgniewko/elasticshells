#include "AverageContactArea.h"

AverageContactArea::AverageContactArea(const char* name, const char* format) : Observer(name, format) {}

AverageContactArea::AverageContactArea(const AverageContactArea& orig) : Observer(orig) {}

AverageContactArea::~AverageContactArea() {}

void AverageContactArea::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double AverageContactArea::observe(const Box& box, std::vector<Shell>& shells, const DomainList& dl)
{
    uint shellsnumber = shells.size();
    double total_contact_a = 0.0;
    double partial_conact_a = 0.0;
    double counter = 0.0;

    for (uint i = 0; i < shellsnumber; i++)
    {
        for (uint j = 0; j < shellsnumber; j++)
        {
            if (i != j)
            {
                partial_conact_a = shells[i].contactArea(shells[j], box);
                total_contact_a += partial_conact_a;

                if (partial_conact_a > 0.0)
                {
                    counter += 1.0;
                }
            }
        }
    }

    if (counter == 0)
    {
        return 0.0;
    }

    return (total_contact_a / counter);
}

DerivedRegister<AverageContactArea> AverageContactArea::reg("AverageContactArea");