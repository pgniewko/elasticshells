#ifndef VOLUMEPRESSURE_H
#define	VOLUMEPRESSURE_H

#include <vector>

#include "simulation/Box.h"
#include "Cell.h"

class VolumePressure {
public:
    VolumePressure();
    VolumePressure(const VolumePressure& orig);
    virtual ~VolumePressure();
    
    static double calcPressure(Box&, vector<Cell>&);
private:

};

#endif	/* VOLUMEPRESSURE_H */

