#include "NullObserver.h"

NullObserver::NullObserver(const char* name, const char* format) : Observer(name, format) {}

NullObserver::NullObserver(const NullObserver& orig) : Observer(orig) {}

NullObserver::~NullObserver() {}

double NullObserver::observe(const Box& box, const std::vector<Shell>& shells)
{
    return 0.0;
}

DerivedRegister<NullObserver> NullObserver::reg("NullObserver");

