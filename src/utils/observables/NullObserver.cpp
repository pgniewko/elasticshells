#include "NullObserver.h"

NullObserver::NullObserver(const char* name, const char* format) : Observer(name, format) {}

NullObserver::NullObserver(const NullObserver& orig) : Observer(orig) {}

NullObserver::~NullObserver() {}

void NullObserver::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double NullObserver::observe(const Box& box, std::vector<Cell>& cells, const DomainList& dl)
{
    return 0.0;
}

DerivedRegister<NullObserver> NullObserver::reg("NullObserver");

