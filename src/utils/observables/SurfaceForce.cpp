#include "SurfaceForce.h"

SurfaceForce::SurfaceForce(const char* name, const char* format) : Observer(name, format) {}

SurfaceForce::SurfaceForce(const SurfaceForce& orig) : Observer(orig) {}

SurfaceForce::~SurfaceForce() {}

void SurfaceForce::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double SurfaceForce::observe(const Box& box, std::vector<Cell>& cells, const DomainList& dl)
{
    return SurfaceForce::calcTotalForce(box, cells);
}

double SurfaceForce::calcTotalForce(const Box& box, std::vector<Cell>& cells)
{
    if (box.pbc)
    {
        return 0.0;
    }

    uint numOfCells = cells.size();
    double totalForce = 0.0;

    for (uint i = 0; i < numOfCells; i++)
    {
        totalForce += cells[i].contactForceSF(box);
    }

    return totalForce;
}

DerivedRegister<SurfaceForce> SurfaceForce::reg("SurfaceForce");
