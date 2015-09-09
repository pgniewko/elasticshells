#include "CellCellStress.h"

CellCellStress::CellCellStress(const char* name, const char* format) : Observer(name, format) {}

CellCellStress::CellCellStress(const CellCellStress& orig) : Observer(orig) {}

CellCellStress::~CellCellStress() {}

void CellCellStress::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double CellCellStress::observe(const Box& box, std::vector<Cell>& cells)
{
    uint cellsnumber = cells.size();
    double contact_force = 0.0;
    double contact_area = 0.0;
    double average_stress = 0.0;
    int counter = 0;

    for (uint i = 0; i < cellsnumber; i++)
    {
        for (uint j = 0; j < cellsnumber; j++)
        {
            if (i != j)
            {
                //contact_force = cells[i].contactForceNew(cells[j], box);
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
        average_stress /= (double)counter;
    }
    else
    {
        average_stress = 0.0;
    }

    return average_stress;
}

DerivedRegister<CellCellStress> CellCellStress::reg("CellCellStress");