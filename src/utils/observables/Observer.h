#ifndef OBSERVER_H
#define	OBSERVER_H

#include <string.h>
#include <stdio.h>
#include <vector>
#include <cstdarg>


#include "Shell.h"
#include "Environment.h"
#include "simulation/Box.h"
#include "simulation/elements.h"
#include "force/ForcesCalculator.h"

class Observer
{
    public:
        explicit Observer(const char*, const char*);
        Observer(const Observer& orig);
        virtual ~Observer();

        virtual void set_params(const int, std::vector<std::string>) = 0;
        virtual double observe(const Box&, const std::vector<Shell>&) = 0;

        const char* getFormat();
        const char* getName();
        
        //const bool is_per_shell() {return per_shell_observer;}
        
        //const int get_i() {return i_param;}
        //const int get_d() {return d_param;}
        
        //void set_i(int new_i) {i_param = new_i;}
        //void set_i(double new_d) {d_param = new_d;}

    protected:

        void create_shells_image(const Box&, const std::vector<Shell>&);

        void copy_shells_data(const Box&, const std::vector<Shell>&);

        bool is_in_contact(int, int);

        const std::string observer_name;
        const std::string output_format;
        int i_param;
        double d_param;

        bool image_not_created = true;
        
        //bool per_shell_observer = false;
        
        std::vector<double> xyz;
        std::vector<double> forces;
        ForcesCalculator fc;

        std::vector<element> elements;
        std::vector<hinge> hinges;
        std::vector<object_map> vs_map;
        std::map<object_map, int> inv_vs_map;

        std::vector<object_map> ts_map;
        std::map<object_map, int> inv_ts_map;

        std::vector<object_map> hs_map;
        std::map<object_map, int> inv_hs_map;

        std::vector<double> turgors;

        std::vector<std::vector<int> > graph;

        pairs_t contacts;
        
        static int MAX_M;
};

template<typename T>
Observer* createT(const char* n, const char* t)
{
    return new T (n, t);
}

struct ObserverFactory
{
        typedef std::map<std::string, Observer * (*)(const char*, const char*)> map_type;

        static Observer* createInstance(std::string const& s, const char* n, const char* t)
        {
            map_type::iterator it = getMap()->find(s);

            if (it == getMap()->end())
            {
                return 0;
            }

            return it->second(n, t);
        }

    protected:
        static map_type* getMap()
        {
            // never delete'ed. (exist until program termination)
            // because we can't guarantee correct destruction order
            if (!map)
            {
                map = new map_type;
            }

            return map;
        }

    public:
        static map_type* map;

};

template<typename T>
struct DerivedRegister : ObserverFactory
{
    DerivedRegister(std::string const& s)
    {
        getMap()->insert( std::make_pair(s, &createT<T>) );
    }
};

#endif	/* OBSERVER_H */