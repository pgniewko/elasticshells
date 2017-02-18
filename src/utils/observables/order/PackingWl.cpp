#include "PackingWl.h"

PackingWl::PackingWl(const char* name, const char* format) : Observer(name, format) 
{}

PackingWl::PackingWl(const PackingWl& orig) : Observer(orig)
{}

PackingWl::~PackingWl() 
{}

void PackingWl::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
    d_param = strtod(args_[ num + 1 ].c_str(), NULL);
}

double PackingWl::observe(const Box& box, std::vector<Cell>& cells)
{
    return 0.0;
}


DerivedRegister<PackingWl> PackingWl::reg("PackingWl");