#ifndef VOLUMEFRACTION_H
#define	VOLUMEFRACTION_H

#include <vector>

#include "simulation/Box.h"
#include "Cell.h"

class VolumeFraction
{
    public:
        VolumeFraction();
        VolumeFraction(const VolumeFraction& orig);
        virtual ~VolumeFraction();
        static double caclVolumeFraction(Box&, std::vector<Cell>&, double = 0.0);
        static double caclCellsVolume(std::vector<Cell>&, double = 0.0);
    private:

};

#endif	/* VOLUMEFRACTION_H */