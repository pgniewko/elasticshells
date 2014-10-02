#ifndef SURFACEPRESSURE_H
#define	SURFACEPRESSURE_H

#include <vector>

#include "simulation/Box.h"
#include "Cell.h"

class SurfacePressure {
public:
    SurfacePressure();
    SurfacePressure(const SurfacePressure& orig);
    virtual ~SurfacePressure();
    
    static double calcPressure(Box&, vector<Cell>&);
private:

};

#endif	/* SURFACEPRESSURE_H */

