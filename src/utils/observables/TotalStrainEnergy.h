#ifndef SURFACESTRAINENERGY_H
#define	SURFACESTRAINENERGY_H

#include "utils/observables/Observer.h"

class TotalStrainEnergy : public Observer
{
    public:
        explicit TotalStrainEnergy(const char*, const char*);
        TotalStrainEnergy(const TotalStrainEnergy& orig);
        virtual ~TotalStrainEnergy();

        double observe(const Box&, std::vector<Cell>&);
        void set_params(const int, std::vector<std::string>);

    private:
        static DerivedRegister<TotalStrainEnergy> reg;

};

#endif	/* SURFACESTRAINENERGY_H */