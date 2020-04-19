#include "PerShellTurgor.h"

PerShellTurgor::PerShellTurgor (const char* name, const char* format) : Observer(name, format) {}

PerShellTurgor::PerShellTurgor(const PerShellTurgor& orig) : Observer(orig) {}

PerShellTurgor::~PerShellTurgor() {}

void PerShellTurgor::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
}

/**
 * Return the i_param-th shell's turgor.
 * 
 * @param box
 * @param shells
 * @return 
 */
double PerShellTurgor::observe(const Box& box, const std::vector<Shell>& shells)
{
    if ((uint)i_param < shells.size())
    {
        return shells[i_param].get_turgor();
    }
    
    return -1;
}

DerivedRegister<PerShellTurgor> PerShellTurgor::reg("PerShellTurgor");

