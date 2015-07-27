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

template<typename T> Observer* createT(const char* n, const char* t) 
{ 
//    std::cout << "inside createT; " << "n=" << n << " t=" << t  << std::endl;
    return new T (n,t);
}

struct ObserverFactory 
{
    typedef std::map<std::string, Observer*(*)(const char*, const char*)> map_type;

    static Observer * createInstance(std::string const& s, const char* n, const char* t) 
    {
        std::cout << "n= " << n << " t=" << t << std::endl;
        map_type::iterator it = getMap()->find(s);
        if(it == getMap()->end())
            return 0;
        
        //std::cout << "it->first=" << it->first << std::endl;
        //std::cout << "it->second=" << it->second << std::endl;
        //(*it->second)("A","B");
        //return (*it->second);//("A","B");
        //std::cout << it->second() << std::endl;
        return it->second(n,t);
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

public:
    static map_type * map;
};

template<typename T>
struct DerivedRegister : ObserverFactory 
{ 
    DerivedRegister(std::string const& s)
    { 
        //std::cout  << "DerivedRegister:jestem tutaj" << std::endl;
        //std::cout  <<  &createT<T> << std::endl;
        //std::cout  <<  &createT<T> << std::endl;
        getMap()->insert( std::make_pair(s, &createT<T>) );
        //std::cout << "dodalem pare" << std::endl;
    }
};

#endif	/* OBSERVER_H */