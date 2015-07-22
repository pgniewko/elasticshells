#ifndef SURFACEFORCE_H
#define	SURFACEFORCE_H

#include <vector>

#include "Cell.h"
#include "simulation/Box.h"
#include "force/HertzianRepulsion.h"


class SurfaceForce
{
    public:
        SurfaceForce();
        SurfaceForce(const SurfaceForce& orig);
        virtual ~SurfaceForce();
        static double calcForces(Box&, std::vector<Cell>&);
    private:

};
#endif	/* SURFACEFORCE_H */
