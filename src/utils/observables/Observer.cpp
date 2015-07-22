#include "Observer.h"

Observer::Observer() {}

Observer::Observer(const Observer& orig) {}

Observer::~Observer() {}

const char* Observer::getFormat()
{
    return output_format;
}

const char* Observer::getName()
{
    return observer_name;
}
