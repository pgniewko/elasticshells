#include "PerShellVolume.h"

PerShellVolume::PerShellVolume (const char* name, const char* format) : Observer(name, format) {}

PerShellVolume::PerShellVolume(const PerShellVolume& orig) : Observer(orig) {}

PerShellVolume::~PerShellVolume() {}

double PerShellVolume::observe(const Box& box, const std::vector<Shell>& shells)
{
    if ((uint)params[0] < shells.size())
    {
        return shells[params[0]].calc_volume(params[0]);
    }
    return -1;
}

DerivedRegister<PerShellVolume> PerShellVolume::reg("PerShellVolume");
