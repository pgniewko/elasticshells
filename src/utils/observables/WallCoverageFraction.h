#ifndef WALLCOVERAGEFRACTION_H
#define	WALLCOVERAGEFRACTION_H

#include <vector>

#include "Cell.h"
#include "simulation/Box.h"

class WallCoverageFraction {
public:
    WallCoverageFraction();
    WallCoverageFraction(const WallCoverageFraction& orig);
    virtual ~WallCoverageFraction();
    static double wallsCoverage(Box&, std::vector<Cell>&);
private:

};

#endif	/* WALLCOVERAGEFRACTION_H */