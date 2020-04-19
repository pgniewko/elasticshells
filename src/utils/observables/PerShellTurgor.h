#ifndef PERSHELLTURGOR_H
#define PERSHELLTURGOR_H

#include "utils/observables/Observer.h"

/**
 * Per shell observer. 
 * It requires only parameters: the shell's index.
 * Configure in the observers.config file as: 
 * PerShellTurgor pst %10.6f 0
 */
class PerShellTurgor : public Observer
{
    public:
        explicit PerShellTurgor(const char*, const char*);
        PerShellTurgor(const PerShellTurgor& orig);
        virtual ~PerShellTurgor();
        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, const std::vector<Shell>&);

    private:
        
        static DerivedRegister<PerShellTurgor> reg;
};

#endif /* PERSHELLTURGOR_H */

