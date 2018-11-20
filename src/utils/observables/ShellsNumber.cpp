#include "ShellsNumber.h"

ShellsNumber::ShellsNumber(const char* name, const char* format) : Observer(name, format) {}

ShellsNumber::ShellsNumber(const ShellsNumber& orig) : Observer(orig) {}

ShellsNumber::~ShellsNumber() {}

void ShellsNumber::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double ShellsNumber::observe(const Box& box, std::vector<Shell>& shells, const DomainList& dl)
{
    return (double) shells.size();
}

DerivedRegister<ShellsNumber> ShellsNumber::reg("ShellsNumber");

