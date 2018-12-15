#include "ShellBoxStress.h"

ShellBoxStress::ShellBoxStress(const char* name, const char* format) : Observer(name, format) {}

ShellBoxStress::ShellBoxStress(const ShellBoxStress& orig) : Observer(orig) {}

ShellBoxStress::~ShellBoxStress() {}

void ShellBoxStress::set_params(const int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
};

double ShellBoxStress::observe(const Box& box, std::vector<Shell>& shells, const DomainList& dl)
{
    if (box.pbc)
    {
        return 0.0;
    }

    uint shellsnumber = shells.size();
    double contact_force = 0.0;
    double contact_area = 0.0;
    double average_stress = 0.0;
    int counter = 0;

    for (uint i = 0; i < shellsnumber; i++)
    {

        contact_force = shells[i].contactForce(box);
        contact_area = shells[i].contactArea(box, d_param);

        if (contact_force > 0 && contact_area > 0)
        {
            average_stress += (contact_force / contact_area);
            counter++;
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

DerivedRegister<ShellBoxStress> ShellBoxStress::reg("ShellBoxStress");