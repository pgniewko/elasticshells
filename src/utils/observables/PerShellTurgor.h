#ifndef PERSHELLTURGOR_H
#define PERSHELLTURGOR_H

#include "utils/observables/Observer.h"

class PerShellTurgor : public Observer
{
    public:
        explicit PerShellTurgor(const char*, const char*);
        PerShellTurgor(const PerShellTurgor& orig);
        virtual ~PerShellTurgor();

        double observe(const Box&, const std::vector<Shell>&);

    private:        
        static DerivedRegister<PerShellTurgor> reg;
};

#endif /* PERSHELLTURGOR_H */

