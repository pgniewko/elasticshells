#ifndef LOGSIMULATION_H
#define	LOGSIMULATION_H

#include <fstream>
#include <vector>

#include "Cell.h"
#include "simulation/Box.h"
#include "utils/observables/SurfacePressure.h"
#include "utils/observables/VolumeFraction.h"
#include "utils/observables/order/QL.h"

class LogSimulation
{
    public:
        LogSimulation(char*);
        LogSimulation(const LogSimulation& orig);
        virtual ~LogSimulation();

        void open();
        void close();
        void dumpState(Box&, std::vector<Cell>&, double, int, int, int);

    private:
        char* logfile;
        FILE* os;

};

#endif	/* LOGSIMULATION_H */

