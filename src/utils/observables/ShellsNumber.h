#ifndef SHELLSNUMBER_H
#define SHELLSNUMBER_H

#include "utils/observables/Observer.h"

class ShellsNumber : public Observer
{
    public:
        explicit ShellsNumber(const char*, const char*);
        ShellsNumber(const ShellsNumber& orig);
        virtual ~ShellsNumber();

        double observe(const Box&, const std::vector<Shell>&);

    private:
        static DerivedRegister<ShellsNumber> reg;
};

#endif /* SHELLSNUMBER_H */

