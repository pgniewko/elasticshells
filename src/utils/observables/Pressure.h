#ifndef PRESSURE_H
#define	PRESSURE_H

#include "utils/observables/Observer.h"
#include "utils/Logger.h"

class Pressure : public Observer
{
    public:
        explicit Pressure(const char*, const char*);
        Pressure(const Pressure& orig);
        virtual ~Pressure();

        double observe(const Box&, const std::vector<Shell>&);

    private:
        double total_force(const Box&, const std::vector<Shell>&);
        static DerivedRegister<Pressure> reg;
        static utils::Logger sp_log;
        bool info_not_printed = true;
};

#endif	/* PRESSURE_H */

