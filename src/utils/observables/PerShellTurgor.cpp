#include "PerShellTurgor.h"

PerShellTurgor::PerShellTurgor (const char* name, const char* format) : Observer(name, format) {}

PerShellTurgor::PerShellTurgor(const PerShellTurgor& orig) : Observer(orig) {}

PerShellTurgor::~PerShellTurgor() {}

double PerShellTurgor::observe(const Box& box, const std::vector<Shell>& shells)
{
    if ((uint)params[0] < shells.size())
    {
        return shells[params[0]].get_turgor();
    }
    return -1;
}

DerivedRegister<PerShellTurgor> PerShellTurgor::reg("PerShellTurgor");

