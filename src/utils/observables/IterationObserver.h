#ifndef ITERATIONOBSERVER_H
#define ITERATIONOBSERVER_H

#include "simulation/integrators/Integrator.h"
#include "utils/observables/Observer.h"

class IterationObserver : public Observer
{
    public:
        explicit IterationObserver(const char*, const char*);
        IterationObserver(const IterationObserver& orig);
        virtual ~IterationObserver();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, const std::vector<Shell>&);

    private:
        static DerivedRegister<IterationObserver> reg;
};

#endif /* ITERATIONOBSERVER_H */

