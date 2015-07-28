#ifndef VOLUMEFRACTION_H
#define	VOLUMEFRACTION_H

#include "utils/observables/Observer.h"

class VolumeFraction : public Observer
{
    public:
        VolumeFraction(const char*, const char*);
        VolumeFraction(const VolumeFraction& orig);
        virtual ~VolumeFraction();

        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);
        double observe(Box&, std::vector<Cell>&);

    private:
        double calcVolumeFraction(Box&, std::vector<Cell>&);
        double calcCellsVolume(std::vector<Cell>&);

        static DerivedRegister<VolumeFraction> reg;
        double rv;
};

#endif	/* VOLUMEFRACTION_H */