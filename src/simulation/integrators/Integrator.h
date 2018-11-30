#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "Environment.h"
#include "simulation/Simulator.h"

class Simulator;

class Integrator
{
    public:
        Integrator();
        explicit Integrator(Simulator*); //, char*);
        Integrator(const Integrator& orig);
        virtual ~Integrator();

        void integrate(Simulator* s);
        void resetParams(Simulator* s);

    private:

        void (Integrator::*integrator)(Simulator*) = 0;

//        void setIntegrator(void (Integrator::*functoall)(Simulator* s));
//        void setIntegrator(Simulator* s, char*);

        void fireIntegrator(Simulator*);
        void _vv(Simulator*);

        static int FIRE_Nmin;
        static int FIRE_N;
        static double FIRE_DT;
        static double FIRE_ALPHA;
        static double FIRE_DTMAX;

        static utils::Logger integrator_logs;

};

#endif /* INTEGRATOR_H */

