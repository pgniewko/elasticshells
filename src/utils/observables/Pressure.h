#ifndef PRESSURE_H
#define	PRESSURE_H

//#include "utils/observables/SurfaceForce.h"
#include "utils/observables/Observer.h"
#include "utils/Logger.h"

class Pressure : public Observer
{
    public:
        explicit Pressure(const char*, const char*);
        Pressure(const Pressure& orig);
        virtual ~Pressure();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, const std::vector<Shell>&);

    private:
        double total_force(const Box&, const std::vector<Shell>&);
        static DerivedRegister<Pressure> reg;
        static utils::Logger sp_log;
};

#endif	/* PRESSURE_H */
