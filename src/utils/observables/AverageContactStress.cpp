#include "AverageContactStress.h"

AverageContactStress::AverageContactStress() {}

AverageContactStress::AverageContactStress(const AverageContactStress& orig) {}

AverageContactStress::~AverageContactStress() {}

double AverageContactStress::caclContactStress(Box& box, std::vector<Cell>& cells)
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
                contact_force = cells[i].contactForce(cells[j], box);
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