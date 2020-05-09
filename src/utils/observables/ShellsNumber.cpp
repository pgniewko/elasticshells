#include "ShellsNumber.h"

ShellsNumber::ShellsNumber(const char* name, const char* format) : Observer(name, format) {}

ShellsNumber::ShellsNumber(const ShellsNumber& orig) : Observer(orig) {}

ShellsNumber::~ShellsNumber() {}

double ShellsNumber::observe(const Box& box, const std::vector<Shell>& shells)
{
    return (double) shells.size();
}

DerivedRegister<ShellsNumber> ShellsNumber::reg("ShellsNumber");

