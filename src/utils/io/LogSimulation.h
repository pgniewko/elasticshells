#ifndef LOGSIMULATION_H
#define	LOGSIMULATION_H

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <cstring>

#include "Cell.h"
#include "utils/Logger.h"
#include "utils/utils.h"
#include "simulation/Box.h"
#include "simulation/DomainList.h"
#include "utils/observables/Observer.h"

class LogSimulation
{
    public:
        LogSimulation(std::string, std::string);
        LogSimulation(const LogSimulation& orig);
        virtual ~LogSimulation();

        void open();
        void close();

        void registerObservers();
        void printHeader();
        void dumpState(Box&, std::vector<Cell>&, const DomainList&);
        
        const std::string getFileName() const;
        
        const std::vector<std::string> readTurgorsFile() const;

    private:
        std::string logfile;
        std::string configfile;
        FILE* os;
        std::vector<Observer*> observers;
        std::vector<std::string> readConfigFile();

        static utils::Logger log_logger;
};

#endif	/* LOGSIMULATION_H */