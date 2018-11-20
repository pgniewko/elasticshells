#ifndef SHELLBOXSTRESS_H
#define	SHELLBOXSTRESS_H

#include "utils/observables/Observer.h"

class ShellBoxStress : public Observer
{
    public:
        explicit ShellBoxStress(const char*, const char*);
        ShellBoxStress(const ShellBoxStress& orig);
        virtual ~ShellBoxStress();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Shell>&, const DomainList&);

    private:
        static DerivedRegister<ShellBoxStress> reg;

};

#endif	/* SHELLBOXSTRESS_H */

