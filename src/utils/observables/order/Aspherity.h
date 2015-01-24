#ifndef ASPHERITY_H
#define	ASPHERITY_H

#include "Environment.h"
#include "Cell.h"
#include "geometry/Vector3D.h"

class Aspherity
{
    public:
        Aspherity();
        Aspherity(const Aspherity& orig);
        virtual ~Aspherity();
        static double calcAspherity(Cell);
    private:

};

#endif	/* ASPHERITY_H */

