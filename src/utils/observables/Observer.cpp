#include "Observer.h"

Observer::Observer(const char* name, const char* format) : 
observer_name(name), output_format(format), i_param(0), d_param(0.0)
{}

Observer::Observer(const Observer& orig) : 
observer_name(orig.observer_name), output_format(orig.output_format), i_param(orig.i_param), d_param(orig.d_param)
{};

Observer::~Observer()
{}

const char* Observer::getFormat()
{
    return output_format.c_str();
}

const char* Observer::getName()
{
    return observer_name.c_str();
}

ObserverFactory::map_type* ObserverFactory::map = NULL;