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
        void reset_params(Simulator* s);
        void set_n(uint n);
        
        static void increment_iter() {ITER_NUM = ITER_NUM + 1;}
        static void increment_total_iter() {TOTAL_ITER_NUM = TOTAL_ITER_NUM + 1;}

        static void reset_iter() {ITER_NUM = 0;}

        static int get_iter_num() {return ITER_NUM;}
        static int get_total_iter_num() {return TOTAL_ITER_NUM;}

        
    private:

        void (Integrator::*integrator)(Simulator*) = 0;

        void fire_integrator(Simulator*);
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
        
        static int ITER_NUM;
        static int TOTAL_ITER_NUM;
};

#endif /* INTEGRATOR_H */

