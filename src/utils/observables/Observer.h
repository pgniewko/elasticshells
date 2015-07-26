#ifndef OBSERVER_H
#define	OBSERVER_H

#include <vector>
#include <cstdarg>

#include "Cell.h"
#include "simulation/Box.h"

class Observer 
{
public:
    Observer(const char*, const char*);
    Observer(const Observer& orig);
    virtual ~Observer();
    virtual double observe(Box&, std::vector<Cell>&) =0;
    virtual void set_params(int, ...) =0;

    const char* getFormat();
    const char* getName();
protected:
    const char* observer_name;
    const char* output_format;

};

template<typename T> Observer * createT(const char* n, const char* t) 
{ 
    return new T (n,t);
}

struct ObserverFactory 
{
    typedef std::map<std::string, Observer*(*)()> map_type;

    static Observer * createInstance(std::string const& s) 
    {
        map_type::iterator it = getMap()->find(s);
        if(it == getMap()->end())
            return 0;
        return it->second();
    }

protected:
    static map_type * getMap() 
    {
        // never delete'ed. (exist until program termination)
        // because we can't guarantee correct destruction order 
        if(!map) 
        {
            map = new map_type;
        }
        return map;
    }

private:
    static map_type * map;
};

template<typename T>
struct DerivedRegister : ObserverFactory 
{ 
    DerivedRegister(std::string const& s)
    { 
        std::cout  << "DerivedRegister:jestem tutaj" << std::endl;
        getMap()->insert(std::make_pair(s, &createT<T>));
    }
};

#endif	/* OBSERVER_H */