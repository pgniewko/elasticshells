#include "IntraCellNB.h"

IntraCellNB::IntraCellNB (const char* name, const char* format) : Observer(name, format) {}

IntraCellNB::IntraCellNB (const IntraCellNB& orig) : Observer(orig) {}

IntraCellNB::~IntraCellNB () {}

void IntraCellNB::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double IntraCellNB::observe(const Box& box, std::vector<Cell>& cells)
{

    double nb_energy = 0.0;
    uint cellsnumber = cells.size();

    for (uint i = 0; i < cellsnumber; i++)
    {
        nb_energy += cells[i].nbIntra(box);
    }

    return nb_energy;
}

DerivedRegister<IntraCellNB> IntraCellNB::reg("IntraCellNB");
