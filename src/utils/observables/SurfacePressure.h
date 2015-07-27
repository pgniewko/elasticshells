#ifndef SURFACEPRESSURE_H
#define	SURFACEPRESSURE_H

#include "Cell.h"
#include "simulation/Box.h"
#include "utils/observables/SurfaceForce.h"
#include "utils/observables/Observer.h"


class SurfacePressure : public Observer
{
    public:
        SurfacePressure(const char*, const char*);
        SurfacePressure(const SurfacePressure& orig);
        virtual ~SurfacePressure();
        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);
        double observe(Box&, std::vector<Cell>&);
        
    private:
        static DerivedRegister<SurfacePressure> reg;
        double rv;
};

#endif	/* SURFACEPRESSURE_H */

