#ifndef VOLUMEPRESSURE_H
#define	VOLUMEPRESSURE_H

#include <vector>

#include "Cell.h"
#include "geometry/Vector3D.h"
#include "simulation/Box.h"

class VolumePressure {
public:
    VolumePressure();
    VolumePressure(const VolumePressure& orig);
    virtual ~VolumePressure();
    
    static double calcPressure(Box&, std::vector<Cell>&, int nbhandler);
private:
    static Vector3D getImage(Box&, const Vector3D&);

};

#endif	/* VOLUMEPRESSURE_H */

