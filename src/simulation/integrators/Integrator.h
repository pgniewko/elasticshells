#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "Environment.h"
#include "simulation/Simulator.h"

class Simulator;

class Integrator
{
    public:
        Integrator();
        explicit Integrator(Simulator*);
        Integrator(const Integrator& orig);
        virtual ~Integrator();

        void integrate(Simulator* s);
        void resetParams(Simulator* s);
        void set_n(uint n);

    private:

        void (Integrator::*integrator)(Simulator*) = 0;

        void fireIntegrator(Simulator*);
        void _vv(Simulator*);
        
        uint n;
        std::vector<double> vel;
        std::vector<double> f_p;

        static int FIRE_Nmin;
        static int FIRE_N;
        static double FIRE_DT;
        static double FIRE_ALPHA;
        static double FIRE_DTMAX;
        static utils::Logger integrator_logs;
};

#endif /* INTEGRATOR_H */

