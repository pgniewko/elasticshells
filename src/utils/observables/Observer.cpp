#include "Observer.h"

Observer::Observer(const char* name, const char* format) : observer_name(name), output_format(format) 
{
//    std::cout << "Observer:Buduja mnie"<<std::endl;
    std::cout << "moje imie: " << observer_name << std::endl;
    std::cout << "moj format: " << output_format << std::endl;
}

Observer::Observer(const Observer& orig) : observer_name(orig.observer_name), output_format(orig.output_format)
{
//    std::cout <<"Obseerver: Kopjuja mnie" << std::endl;
}

Observer::~Observer() 
{
//    std::cout <<"Obseerver: Niszcza mnie" << std::endl;
}

const char* Observer::getFormat()
{
    return output_format;
}

const char* Observer::getName()
{
    return observer_name;
}

ObserverFactory::map_type* ObserverFactory::map = NULL;