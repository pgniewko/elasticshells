#include "LogSimulation.h"

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) 
    {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    
    return elems;
}


LogSimulation::LogSimulation(char* lf, char* cf)
{
    logfile = lf;
    configfile = cf;
}

LogSimulation::LogSimulation(const LogSimulation& orig) : logfile(orig.logfile), configfile(orig.configfile)
{
}

LogSimulation::~LogSimulation() 
{
    for (uint i = 0; i < observers.size(); i++)
    {
        delete observers[i];
    }
}

void LogSimulation::open()
{
    os = fopen(logfile, "w");
}

void LogSimulation::close()
{
    fclose(os);
}

std::vector<std::string> LogSimulation::readConfigFile()
{
    std::ifstream cfile;
    cfile.open(configfile);
    
    std::vector<std::string> list;
    std::string line;
    
    
    while( std::getline (cfile, line) )
    {
        if ( !(line.at(0) == '#') && !(line.at(0) == ' ') )
        {
            list.push_back(line);
        }
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
                //std::cout << single_line[0] <<std::endl;
                obs_obj->set_params(3, single_line);
            }
            
            std::cout << obs_obj->getName()<< " " << obs_obj->getFormat() << std::endl;
            observers.push_back( obs_obj );
        }
    }
}

void LogSimulation::printHeader()
{
    fprintf(os,"#");
    
    for(std::vector<Observer*>::iterator it = observers.begin(); it != observers.end(); ++it)
    {
        fprintf(os, std::strcat ( (*it)->getName(), " ") );
    }
    
   fprintf(os,"\n");
   //fflush(os);
}

void LogSimulation::dumpState(Box& box, std::vector<Cell>& cells)
{
    //std::cout << "dumpState" << std::endl;
    for(std::vector<Observer*>::iterator it = observers.begin(); it != observers.end(); ++it)
    {
       std::cout << (*it)->getName() << " "<< (*it)->getFormat()<< "observe="<< (*it)->observe(box, cells) << std::endl;
       
       fprintf(os, (*it)->getFormat(), (*it)->observe(box, cells) );
       fprintf(os," ");
    }
    fprintf(os, "\n");
    //fflush(os);
}