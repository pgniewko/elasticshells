#ifndef LOGSIMULATION_H
#define	LOGSIMULATION_H

#include <fstream>
#include <vector>
#include <string>
#include <memory>

#include "Cell.h"
#include "simulation/Box.h"
#include "utils/observables/SurfacePressure.h"
#include "utils/observables/SurfaceForce.h"
#include "utils/observables/VolumeFraction.h"
#include "utils/observables/order/QL.h"
#include "utils/observables/AverageContactStress.h"
#include "utils/observables/WallCoverageFraction.h"
#include "utils/observables/TotalCellsArea.h"
#include "utils/observables/SurfaceStrainEnergy.h"
#include "utils/observables/AverageTurgor.h"
#include "utils/observables/Observer.h"

class LogSimulation
{
    public:
        LogSimulation(char*, char*);
        LogSimulation(const LogSimulation& orig);
        virtual ~LogSimulation();

        void open();
        void close();
        
        void registerObservers();
        void printHeader();
        void dumpState(Box&, std::vector<Cell>&, double, int, int, int);

    private:
        char* logfile;
        char* configfile;
        FILE* os;
        std::vector<Observer*> observers;
        
        void readConfigFile();

};

#endif	/* LOGSIMULATION_H */

