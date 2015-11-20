#include "AverageActiveArea.h"

AverageActiveArea::AverageActiveArea(const char* name, const char* format) : Observer(name, format) {}

AverageActiveArea::AverageActiveArea(const AverageActiveArea& orig) : Observer(orig) {}

AverageActiveArea::~AverageActiveArea() {}

void AverageActiveArea::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double AverageActiveArea::observe(const Box& box, std::vector<Cell>& cells)
{
    uint cellsnumber = cells.size();
    double total_active_a = 0.0;
    double counter = 0.0;
    
    for (uint i = 0; i < cellsnumber; i++)
    {
        total_active_a += cells[i].activeArea(box, cells);
        counter += 1.0;
    }
    
    return (total_active_a /= counter);
}

DerivedRegister<AverageActiveArea> AverageActiveArea::reg("AverageActiveArea");