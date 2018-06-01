#include "BoxSizeObserver.h"

BoxSizeObserver::BoxSizeObserver(const char* name, const char* format) : Observer(name, format) {}

BoxSizeObserver::BoxSizeObserver(const BoxSizeObserver& orig) : Observer(orig) {}

BoxSizeObserver::~BoxSizeObserver() {}

void BoxSizeObserver::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
};

double BoxSizeObserver::observe(const Box& box, std::vector<Shell>& cells, const DomainList& dl)
{
    switch (i_param)
    {
        case 1:
            return box.getX();

        //break;

        case 2:
            return box.getY();

        //break;

        case 3:
            return box.getZ();

        //break;

        default:
            return 0.0;
            //break;
    }

    return 0.0;
}

DerivedRegister<BoxSizeObserver> BoxSizeObserver::reg("BoxSizeObserver");