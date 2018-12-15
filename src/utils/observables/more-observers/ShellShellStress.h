#ifndef AVERAGECONTACTSTRESSNEW_H
#define	AVERAGECONTACTSTRESSNEW_H

#include "utils/observables/Observer.h"

class ShellShellStress : public Observer
{
    public:
        explicit ShellShellStress(const char*, const char*);
        ShellShellStress(const ShellShellStress& orig);
        virtual ~ShellShellStress();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Shell>&, const DomainList&);

    private:
        static DerivedRegister<ShellShellStress> reg;

};

#endif	/* AVERAGECONTACTSTRESSNEW_H */

