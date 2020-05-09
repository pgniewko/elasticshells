#ifndef PERSHELLVOLUME_H
#define PERSHELLVOLUME_H

#include "utils/observables/Observer.h"

class PerShellVolume : public Observer
{
    public:
        explicit PerShellVolume(const char*, const char*);
        PerShellVolume(const PerShellVolume& orig);
        virtual ~PerShellVolume();

        double observe(const Box&, const std::vector<Shell>&);

    private: 
        static DerivedRegister<PerShellVolume> reg;
};

#endif /* PERSHELLVOLUME_H */

