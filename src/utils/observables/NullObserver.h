#ifndef NULLOBSERVER_H
#define NULLOBSERVER_H

#include "utils/observables/Observer.h"

class NullObserver : public Observer
{
    public:
        NullObserver(const char*, const char*);
        NullObserver(const NullObserver& orig);
        virtual ~NullObserver();
        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, const std::vector<Shell>&);

    private:
        static DerivedRegister<NullObserver> reg;

};

#endif /* NULLOBSERVER_H */

