#ifndef LOGSIMULATION_H
#define	LOGSIMULATION_H

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <cstring>

#include "Shell.h"
#include "utils/Logger.h"
#include "utils/utils.h"
#include "simulation/Box.h"
#include "utils/observables/Observer.h"

class LogSimulation
{
    public:
        LogSimulation(std::string, std::string);
        LogSimulation(const LogSimulation& orig);
        virtual ~LogSimulation();

        void open();
        void close();

        void register_observers();
        void print_header();
        void dump_state(const Box&, const std::vector<Shell>&);

        const std::string get_file_name() const;

        const std::vector<std::string> read_turgors_file() const;

    private:
        std::string logfile;
        std::string configfile;
        FILE* os;
        std::vector<Observer*> observers;
        std::vector<std::string> read_configuration_file();

        static utils::Logger log_logger;
};

#endif	/* LOGSIMULATION_H */