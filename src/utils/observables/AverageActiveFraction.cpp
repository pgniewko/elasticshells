#include "AverageActiveFraction.h"

AverageActiveFraction::AverageActiveFraction(const char* name, const char* format) : Observer(name, format) {}

AverageActiveFraction::AverageActiveFraction(const AverageActiveFraction& orig) : Observer(orig) {}

AverageActiveFraction::~AverageActiveFraction() {}

void AverageActiveFraction::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
    return;
};

double AverageActiveFraction::observe(const Box& box, std::vector<Cell>& cells)
{
    bool flag = false;
    if (i_param > 0)
        flag = true;
    
    uint cellsnumber = cells.size();
    double total_active_f = 0.0;
    double counter = 0.0;
    
    for (uint i = 0; i < cellsnumber; i++)
    {
        total_active_f += cells[i].activeAreaFraction(box, cells, counter, flag);
        //counter += 1.0;
    }
    
    return (total_active_f /= counter);
}

DerivedRegister<AverageActiveFraction> AverageActiveFraction::reg("AverageActiveFraction");
