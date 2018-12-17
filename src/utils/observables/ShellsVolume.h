#ifndef SHELLSVOLUME_H
#define	SHELLSVOLUME_H

#include <algorithm>
#include "utils/observables/Observer.h"

class ShellsVolume : public Observer
{
    public:
        explicit ShellsVolume(const char*, const char*);
        ShellsVolume(const ShellsVolume& orig);
        virtual ~ShellsVolume();
        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, const std::vector<Shell>&);

    private:
        double intersection(const Shell&, const Shell&);

        static DerivedRegister<ShellsVolume> reg;
};

#endif	/* SHELLSVOLUME_H */

