#include "PerShellVolume.h"

PerShellVolume::PerShellVolume (const char* name, const char* format) : Observer(name, format) {}

PerShellVolume::PerShellVolume(const PerShellVolume& orig) : Observer(orig) {}

PerShellVolume::~PerShellVolume() {}

void PerShellVolume::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
    d_param = strtod(args_[ num + 1 ].c_str(), NULL);
}

/**
 * Return the i_param-th shell's corrected volume.
 * 
 * @param box
 * @param shells
 * @return 
 */
double PerShellVolume::observe(const Box& box, const std::vector<Shell>& shells)
{
    if ((uint)i_param < shells.size())
    {
        return shells[i_param].calc_volume(d_param);
    }
    
    return -1;
}

DerivedRegister<PerShellVolume> PerShellVolume::reg("PerShellVolume");
