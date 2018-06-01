#include "CellBoxStress.h"

CellBoxStress::CellBoxStress(const char* name, const char* format) : Observer(name, format) {}

CellBoxStress::CellBoxStress(const CellBoxStress& orig) : Observer(orig) {}

CellBoxStress::~CellBoxStress() {}

void CellBoxStress::set_params(const int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
};

double CellBoxStress::observe(const Box& box, std::vector<Shell>& cells, const DomainList& dl)
{
    if (box.pbc)
    {
        return 0.0;
    }

    uint cellsnumber = cells.size();
    double contact_force = 0.0;
    double contact_area = 0.0;
    double average_stress = 0.0;
    int counter = 0;

    for (uint i = 0; i < cellsnumber; i++)
    {

        contact_force = cells[i].contactForce(box);
        contact_area = cells[i].contactArea(box, d_param);

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

DerivedRegister<CellBoxStress> CellBoxStress::reg("CellBoxStress");