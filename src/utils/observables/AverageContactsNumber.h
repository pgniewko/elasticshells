#ifndef AVERAGECONTACTSNUMBER_H
#define AVERAGECONTACTSNUMBER_H

#include "utils/observables/Observer.h"

class AverageContactsNumber : public Observer
{
    public:
        explicit AverageContactsNumber(const char*, const char*);
        AverageContactsNumber(const AverageContactsNumber& orig);
        virtual ~AverageContactsNumber();
    
        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&, const DomainList&);
        
    private:
        static DerivedRegister<AverageContactsNumber> reg;
};

#endif /* CONTACTSNUMBER_H */

