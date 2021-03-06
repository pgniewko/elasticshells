#ifndef AVERAGECONTACTSNUMBER_H
#define AVERAGECONTACTSNUMBER_H

#include "utils/observables/Observer.h"

class AverageContactsNumber : public Observer
{
    public:
        explicit AverageContactsNumber(const char*, const char*);
        AverageContactsNumber(const AverageContactsNumber& orig);
        virtual ~AverageContactsNumber();

        double observe(const Box&, const std::vector<Shell>&);

    private:
        static DerivedRegister<AverageContactsNumber> reg;
};

#endif /* CONTACTSNUMBER_H */

