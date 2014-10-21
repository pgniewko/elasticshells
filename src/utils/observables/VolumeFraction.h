#ifndef VOLUMEFRACTION_H
#define	VOLUMEFRACTION_H

#include <vector>

#include "simulation/Box.h"
#include "Cell.h"

class VolumeFraction {
public:
    VolumeFraction();
    VolumeFraction(const VolumeFraction& orig);
    virtual ~VolumeFraction();
    static double caclVolumeFraction(Box&, vector<Cell>&);
    static double caclCellsVolume(vector<Cell>&);
private:

};

#endif	/* VOLUMEFRACTION_H */

