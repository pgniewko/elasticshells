#include "IterationObserver.h"

IterationObserver::IterationObserver(const char* name, const char* format) : Observer(name, format) {}

IterationObserver::IterationObserver(const IterationObserver& orig) : Observer(orig) {}

IterationObserver::~IterationObserver() {}

void IterationObserver::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double IterationObserver::observe(const Box& box, const std::vector<Shell>& shells)
{
    return (double) Integrator::get_total_iter_num();
}

DerivedRegister<IterationObserver> IterationObserver::reg("IterationObserver");

