#include "LogSimulation.h"

utils::Logger LogSimulation::log_logger("log_logger");

LogSimulation::LogSimulation(std::string lf, std::string cf) : logfile(lf), configfile(cf) {}

LogSimulation::LogSimulation(const LogSimulation& orig) : logfile(orig.logfile), configfile(orig.configfile) {}

LogSimulation::~LogSimulation()
{
    for (uint i = 0; i < observers.size(); i++)
    {
        delete observers[i];
    }
}

void LogSimulation::open()
{
    os = fopen(logfile.c_str(), "a");

    if ( os == NULL )
    {
        os = fopen(logfile.c_str(), "w");
    }

    if ( os == NULL )
    {
        log_logger << utils::LogLevel::WARNING << "Can not open file:<<" << logfile << "for writing.\n";
    }

    return;
}

void LogSimulation::close()
{
    if ( os != NULL )
    {
        fclose(os);
    }
}

std::vector<std::string> LogSimulation::readConfigFile()
{
    std::ifstream cfile;
    cfile.open(configfile, std::ifstream::in);
    std::vector<std::string> list;
    std::string line;

    if ( cfile.is_open() )
    {
        while ( std::getline (cfile, line) )
        {
            if ( !(line.at(0) == '#') && !(line.at(0) == ' ') )
            {
                list.push_back(line);
            }
        }
    }
    else
    {
        log_logger << utils::LogLevel::WARNING << "Observers configuration file COULD NOT BE FOUND" << "\n";
    }

    cfile.close();
    return list;
}

void LogSimulation::registerObservers()
{
    std::vector<std::string> list = readConfigFile();
    std::vector<std::string> single_line;

    for (std::vector<std::string>::iterator it = list.begin(); it != list.end(); ++it)
    {
        single_line = split( *it, ' ');

        if (single_line.size() >= 3)
        {
            Observer* obs_obj = ObserverFactory::createInstance( single_line[0], single_line[1].c_str(), single_line[2].c_str() );

            if (single_line.size() > 3)
            {
                obs_obj->set_params(3, single_line);
            }

            observers.push_back( obs_obj );
        }
    }
}

void LogSimulation::printHeader()
{
    fprintf(os, "%s ", "#");

    for (std::vector<Observer*>::iterator it = observers.begin(); it != observers.end(); ++it)
    {
        fprintf(os, "%s ", (*it)->getName());
    }

    fprintf(os, "%s", "\n");
    fflush(os);
}

void LogSimulation::dumpState(Box& box, std::vector<Cell>& cells)
{
    for (std::vector<Observer*>::iterator it = observers.begin(); it != observers.end(); ++it)
    {
        fprintf(os, (*it)->getFormat(), (*it)->observe(box, cells) );
        fprintf(os, "%s", " ");
    }

    fprintf(os, "%s" , "\n");
    fflush(os);
}