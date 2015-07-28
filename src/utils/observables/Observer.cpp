#include "Observer.h"

Observer::Observer(const char* name, const char* format) : observer_name(name), output_format(format) 
{}

Observer::Observer(const Observer& orig) : observer_name(orig.observer_name), output_format(orig.output_format)
{}

Observer::~Observer() 
{}

const char* Observer::getFormat()
{
    return output_format;
}

const char* Observer::getName()
{
    return observer_name;
}

ObserverFactory::map_type* ObserverFactory::map = NULL;