#ifndef SURFACEPRESSURE_H
#define	SURFACEPRESSURE_H

#include <vector>

#include "Cell.h"
#include "simulation/Box.h"
#include "force/HertzianRepulsion.h"

class SurfacePressure {
public:
    SurfacePressure();
    SurfacePressure(const SurfacePressure& orig);
    virtual ~SurfacePressure();
    
    static double calcPressure(Box&, std::vector<Cell>&);
private:

};

#endif	/* SURFACEPRESSURE_H */

