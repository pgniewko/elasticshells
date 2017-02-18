#include "PackingQl.h"

PackingQl::PackingQl(const char* name, const char* format) : Observer(name, format)
{}

PackingQl::PackingQl(const PackingQl& orig) : Observer(orig)
{}

PackingQl::~PackingQl() {
}

void PackingQl::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
    d_param = strtod(args_[ num + 1 ].c_str(), NULL);
}

double PackingQl::observe(const Box& box, std::vector<Cell>& cells)
{
    return 0.0;
}


DerivedRegister<PackingQl> PackingQl::reg("PackingQl");
