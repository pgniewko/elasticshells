#ifndef SURFACEPRESSURE_H
#define	SURFACEPRESSURE_H

#include "utils/observables/SurfaceForce.h"
#include "utils/observables/Observer.h"

class SurfacePressure : public Observer
{
    public:
        SurfacePressure(const char*, const char*);
        SurfacePressure(const SurfacePressure& orig);
        virtual ~SurfacePressure();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);

    private:
        static DerivedRegister<SurfacePressure> reg;
};

#endif	/* SURFACEPRESSURE_H */

