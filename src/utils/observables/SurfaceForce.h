#ifndef SURFACEFORCE_H
#define	SURFACEFORCE_H

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

        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);
        double observe(Box&, std::vector<Cell>&);
        
    private:
        static double calcTotalForce(Box&, std::vector<Cell>&);
        static DerivedRegister<SurfaceForce> reg;

};
#endif	/* SURFACEFORCE_H */
