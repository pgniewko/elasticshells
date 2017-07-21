#ifndef AVERAGECONTACTFORCE_H
#define AVERAGECONTACTFORCE_H

#include "utils/observables/Observer.h"

class AverageContactForce : public Observer
{
    public:
        AverageContactForce(const char*, const char*);
        AverageContactForce(const AverageContactForce& orig);
        virtual ~AverageContactForce();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&, const DomainList&);

    private:
        static DerivedRegister<AverageContactForce> reg;
};

#endif /* AVERAGECONTACTFORCE_H */

