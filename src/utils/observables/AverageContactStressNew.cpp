#include "AverageContactStressNew.h"

AverageContactStressNew::AverageContactStressNew(const char* name, const char* format) : Observer(name, format)
{}

AverageContactStressNew::AverageContactStressNew(const AverageContactStressNew& orig) : Observer(orig)
{}

AverageContactStressNew::~AverageContactStressNew() 
{}

void AverageContactStressNew::set_params(int num, ...)
{
    return;
};

void AverageContactStressNew::set_params(int num, std::vector<std::string> args_)
{
    return;
};

double AverageContactStressNew::observe(Box& box, std::vector<Cell>& cells)
{
    int cellsnumber = cells.size();
    double contact_force = 0.0;
    double contact_area = 0.0;
    double average_stress = 0.0;
    int counter = 0;

    for (int i = 0; i < cellsnumber; i++)
    {
        for (int j = 0; j < cellsnumber; j++)
        {
            if (i != j)
            {
                contact_force = cells[i].contactForceNew(cells[j], box);
                contact_area = cells[i].contactArea(cells[j], box);

                if (contact_force > 0 && contact_area > 0)
                {
                    average_stress += (contact_force / contact_area);
                    counter++;
                }
            }
        }
    }

    if (counter > 0)
    {
        average_stress /= counter;
    }
    else
    {
        average_stress = 0.0;
    }

    return average_stress;
}

DerivedRegister<AverageContactStressNew> AverageContactStressNew::reg("AverageContactStressNew");

