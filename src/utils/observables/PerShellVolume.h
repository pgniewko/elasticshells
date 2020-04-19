#ifndef PERSHELLVOLUME_H
#define PERSHELLVOLUME_H

#include "utils/observables/Observer.h"

/**
 * Per shell volume observer. 
 * It requires two parameters: the shell's index and the volume correction.
 * Configure in the observers.config file as: 
 * PerShellVolume psv %10.6f 0 0.1
 */
class PerShellVolume : public Observer
{
    public:
        explicit PerShellVolume(const char*, const char*);
        PerShellVolume(const PerShellVolume& orig);
        virtual ~PerShellVolume();
        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, const std::vector<Shell>&);

    private:
        
        static DerivedRegister<PerShellVolume> reg;
};

#endif /* PERSHELLVOLUME_H */

