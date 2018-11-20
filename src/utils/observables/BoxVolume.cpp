#include "BoxVolume.h"

BoxVolume::BoxVolume (const char* name, const char* format) : Observer(name, format) {}

BoxVolume::BoxVolume(const BoxVolume& orig) : Observer(orig) {}

BoxVolume::~BoxVolume() {}

void BoxVolume::set_params(const int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
}

double BoxVolume::observe(const Box& box, std::vector<Shell>& shells, const DomainList& dl)
{
    return box.getVolume(d_param);
}

DerivedRegister<BoxVolume> BoxVolume::reg("BoxVolume");

