#ifndef SURFACESTRAINENERGY_H
#define	SURFACESTRAINENERGY_H

#include "utils/observables/Observer.h"

class SurfaceStrainEnergy : public Observer
{
    public:
        SurfaceStrainEnergy(const char*, const char*);
        SurfaceStrainEnergy(const SurfaceStrainEnergy& orig);
        virtual ~SurfaceStrainEnergy();

        double observe(Box&, std::vector<Cell>&);
        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);

    private:
        static DerivedRegister<SurfaceStrainEnergy> reg;

};

#endif	/* SURFACESTRAINENERGY_H */