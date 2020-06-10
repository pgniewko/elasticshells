#include "BoxVolume.h"

BoxVolume::BoxVolume (const char* name, const char* format) : Observer(name, format) {}

BoxVolume::BoxVolume(const BoxVolume& orig) : Observer(orig) {}

BoxVolume::~BoxVolume() {}

double BoxVolume::observe(const Box& box, const std::vector<Shell>& shells)
{
    return box.get_volume(params[0]);
}

DerivedRegister<BoxVolume> BoxVolume::reg("BoxVolume");