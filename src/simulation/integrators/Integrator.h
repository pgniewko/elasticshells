#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "../Simulator.h"

class Integrator 
{
    public:
        explicit Integrator();
        Integrator(const Integrator& orig);
        virtual ~Integrator();
    
        virtual void integrate(Simulator*) = 0;
    private:
        
        Simulator* my_simulator;

};

#endif /* INTEGRATOR_H */

