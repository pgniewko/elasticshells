#ifndef SURFACEFORCE_H
#define	SURFACEFORCE_H

#include "force/HertzianRepulsion.h"
#include "utils/observables/Observer.h"

class Pressure;

class SurfaceForce : public Observer
{
        friend class Pressure;
    public:
        explicit SurfaceForce(const char*, const char*);
        SurfaceForce(const SurfaceForce& orig);
        virtual ~SurfaceForce();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);

    private:
        static double calcTotalForce(const Box&, std::vector<Cell>&);
        static DerivedRegister<SurfaceForce> reg;

};
#endif	/* SURFACEFORCE_H */
