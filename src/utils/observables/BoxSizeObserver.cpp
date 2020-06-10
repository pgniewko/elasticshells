#include "BoxSizeObserver.h"

BoxSizeObserver::BoxSizeObserver(const char* name, const char* format) : Observer(name, format) {}

BoxSizeObserver::BoxSizeObserver(const BoxSizeObserver& orig) : Observer(orig) {}

BoxSizeObserver::~BoxSizeObserver() {}

double BoxSizeObserver::observe(const Box& box, const std::vector<Shell>& shells)
{
    switch ((int)params[0])
    {
        case 1:
            return box.get_x();

        case 2:
            return box.get_y();


        case 3:
            return box.get_z();


        default:
            return 0.0;
    }

    return 0.0;
}

DerivedRegister<BoxSizeObserver> BoxSizeObserver::reg("BoxSizeObserver");