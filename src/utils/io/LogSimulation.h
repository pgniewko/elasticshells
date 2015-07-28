#ifndef LOGSIMULATION_H
#define	LOGSIMULATION_H

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <cstring>

#include "Cell.h"
#include "simulation/Box.h"

#include "utils/observables/Observer.h"

//#include "utils/observables/order/Aspherity.h"
//#include "utils/observables/order/QL.h"
//#include "utils/observables/order/WL.h"

//#include "utils/observables/SurfacePressure.h"
//#include "utils/observables/SurfaceForce.h"
//#include "utils/observables/VolumeFraction.h"
//#include "utils/observables/order/QL.h"
//#include "utils/observables/AverageContactStress.h"
//#include "utils/observables/WallCoverageFraction.h"
//#include "utils/observables/TotalCellsArea.h"
//#include "utils/observables/SurfaceStrainEnergy.h"
//#include "utils/observables/AverageTurgor.h"

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
        void dumpState(Box&, std::vector<Cell>&);

    private:
        char* logfile;
        char* configfile;
        FILE* os;
        std::vector<Observer*> observers;
        std::vector<std::string> readConfigFile();
};

#endif	/* LOGSIMULATION_H */

