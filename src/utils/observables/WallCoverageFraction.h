#ifndef WALLCOVERAGEFRACTION_H
#define	WALLCOVERAGEFRACTION_H

#include "utils/observables/Observer.h"

class WallCoverageFraction : public Observer
{
    public:
        WallCoverageFraction(const char*, const char*);
        WallCoverageFraction(const WallCoverageFraction& orig);
        virtual ~WallCoverageFraction();

        double observe(Box&, std::vector<Cell>&);
        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);

    private:
        static DerivedRegister<WallCoverageFraction> reg;
};

#endif	/* WALLCOVERAGEFRACTION_H */