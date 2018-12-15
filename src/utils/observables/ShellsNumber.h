#ifndef SHELLSNUMBER_H
#define SHELLSNUMBER_H

#include "utils/observables/Observer.h"

class ShellsNumber : public Observer
{
    public:
        explicit ShellsNumber(const char*, const char*);
        ShellsNumber(const ShellsNumber& orig);
        virtual ~ShellsNumber();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, const std::vector<Shell>&);

    private:
        static DerivedRegister<ShellsNumber> reg;

};

#endif /* SHELLSNUMBER_H */

