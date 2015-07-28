#ifndef SURFACEFORCE_H
#define	SURFACEFORCE_H

#include <vector>

#include "force/HertzianRepulsion.h"

#include "utils/observables/Observer.h"

class SurfacePressure;

class SurfaceForce : public Observer
{
        friend class SurfacePressure;
    public:
        SurfaceForce(const char*, const char*);
        SurfaceForce(const SurfaceForce& orig);
        virtual ~SurfaceForce();

        double observe(Box&, std::vector<Cell>&);

        void set_params(int, ...)
        {
            return;
        };
        void set_params(int, std::vector<std::string>)
        {
            return;
        };
    private:
        static double calcTotalForce(Box&, std::vector<Cell>&);
        static DerivedRegister<SurfaceForce> reg;

};
#endif	/* SURFACEFORCE_H */
