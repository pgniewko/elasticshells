#ifndef AVERAGECONTACTSTRESS_H
#define	AVERAGECONTACTSTRESS_H

#include <vector>

#include "Cell.h"
#include "simulation/Box.h"

class AverageContactStress
{
    public:
        AverageContactStress();
        AverageContactStress(const AverageContactStress& orig);
        virtual ~AverageContactStress();
        static double caclContactStress(Box&, std::vector<Cell>&);
    private:

};

#endif	/* AVERAGECONTACTSTRESS_H */