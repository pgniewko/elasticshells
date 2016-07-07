#include "AverageActiveArea.h"

AverageActiveArea::AverageActiveArea(const char* name, const char* format) : Observer(name, format) {}

AverageActiveArea::AverageActiveArea(const AverageActiveArea& orig) : Observer(orig) {}

AverageActiveArea::~AverageActiveArea() {}

void AverageActiveArea::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
    return;
};

double AverageActiveArea::observe(const Box& box, std::vector<Cell>& cells)
{
    bool flag = false;

    if (i_param > 0)
    {
        flag = true;
    }

    uint cellsnumber = cells.size();
    double total_active_a = 0.0;
    double counter = 0.0;

    for (uint i = 0; i < cellsnumber; i++)
    {
        total_active_a += cells[i].activeArea(box, cells, counter, flag);
        //counter += 1.0;
    }

    if (counter == 0.0 )
    {
        return 0.0;
    }

    return (total_active_a /= counter);
}

DerivedRegister<AverageActiveArea> AverageActiveArea::reg("AverageActiveArea");