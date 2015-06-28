#ifndef SURFACESTRAINENERGY_H
#define	SURFACESTRAINENERGY_H

#include <vector>
#include "Cell.h"

class SurfaceStrainEnergy
{
    public:
        SurfaceStrainEnergy();
        SurfaceStrainEnergy(const SurfaceStrainEnergy& orig);
        virtual ~SurfaceStrainEnergy();
        static double calcSurfaceEnergy(std::vector<Cell>&);
    private:

};

#endif	/* SURFACESTRAINENERGY_H */