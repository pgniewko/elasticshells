#include "ShapeIndex.h"

ShapeIndex::ShapeIndex(const char* name, const char* format) : Observer(name, format) {}

ShapeIndex::ShapeIndex(const ShapeIndex& orig) : Observer(orig) {}

ShapeIndex::~ShapeIndex() {}

void ShapeIndex::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double ShapeIndex::observe(const Box& box, std::vector<Shell>& shells, const DomainList& dl)
{
    double q = 0.0;

    double a, v;

    for (uint i = 0; i < shells.size(); i++)
    {
        a = shells[i].calcSurfaceArea(0.0);
        v = shells[i].calcVolume(0.0);;
        q += ( a / pow(v, 2.0 / 3.0) ) ;
    }

    q /= shells.size();
    return q;
}

DerivedRegister<ShapeIndex> ShapeIndex::reg("ShapeIndex");