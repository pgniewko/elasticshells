#ifndef VOLUMEFRACTION_H
#define	VOLUMEFRACTION_H

#include "utils/observables/Observer.h"

class VolumeFraction : public Observer
{
    public:
        VolumeFraction(const char*, const char*);
        VolumeFraction(const VolumeFraction& orig);
        virtual ~VolumeFraction();

        double observe(const Box&, std::vector<Cell>&);
        void set_params(const int, std::vector<std::string>);

    private:
        double calcVolumeFraction(const Box&, std::vector<Cell>&);
        double calcCellsVolume(std::vector<Cell>&);

        static DerivedRegister<VolumeFraction> reg;
};

#endif	/* VOLUMEFRACTION_H */