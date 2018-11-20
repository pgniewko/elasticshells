#ifndef SHELLSVOLUME_H
#define	SHELLSVOLUME_H

#include "utils/observables/Observer.h"
#include <algorithm>    // std::max

class ShellsVolume : public Observer
{
    public:
        explicit ShellsVolume(const char*, const char*);
        ShellsVolume(const ShellsVolume& orig);
        virtual ~ShellsVolume();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Shell>&, const DomainList&);

    private:
        double calcVolumeFraction(const Box&, std::vector<Shell>&);
        double calcShellsVolume(std::vector<Shell>&);

        double sphereSphereIntersection(const Shell&, const Shell&);

        static DerivedRegister<ShellsVolume> reg;
};

#endif	/* SHELLSVOLUME_H */

