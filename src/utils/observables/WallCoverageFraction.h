#ifndef WALLCOVERAGEFRACTION_H
#define	WALLCOVERAGEFRACTION_H

#include "utils/observables/Observer.h"

class WallCoverageFraction : public Observer
{
    public:
        explicit WallCoverageFraction(const char*, const char*);
        WallCoverageFraction(const WallCoverageFraction& orig);
        virtual ~WallCoverageFraction();

        double observe(const Box&, std::vector<Cell>&);
        void set_params(const int, std::vector<std::string>);

    private:
        static DerivedRegister<WallCoverageFraction> reg;
};

#endif	/* WALLCOVERAGEFRACTION_H */