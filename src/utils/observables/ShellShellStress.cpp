#include "ShellShellStress.h"

ShellShellStress::ShellShellStress(const char* name, const char* format) : Observer(name, format) {}

ShellShellStress::ShellShellStress(const ShellShellStress& orig) : Observer(orig) {}

ShellShellStress::~ShellShellStress() {}

void ShellShellStress::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double ShellShellStress::observe(const Box& box, std::vector<Shell>& shells, const DomainList& dl)
{
    uint shellsnumber = shells.size();
    double contact_force = 0.0;
    double contact_area = 0.0;
    double average_stress = 0.0;
    int counter = 0;

    for (uint i = 0; i < shellsnumber; i++)
    {
        for (uint j = 0; j < shellsnumber; j++)
        {
            if (i != j)
            {
                contact_force = shells[i].contactForce(shells[j], box, false);
                contact_area = shells[i].contactArea(shells[j], box);

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
        average_stress /= (double)counter;
    }
    else
    {
        average_stress = 0.0;
    }

    return average_stress;
}

DerivedRegister<ShellShellStress> ShellShellStress::reg("ShellShellStress");