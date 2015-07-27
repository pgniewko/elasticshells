#include "LogSimulation.h"

LogSimulation::LogSimulation(char* lf, char* cf)
{
    logfile = lf;
    configfile = cf;
}

LogSimulation::LogSimulation(const LogSimulation& orig) : logfile(orig.logfile), configfile(orig.configfile) {std::cout << "LogSimulation: Kopjuja mnie"<<std::endl;}

LogSimulation::~LogSimulation() 
{
    std::cout << "LogSimulation: Niszcza mnie"<<std::endl;
    for (uint i = 0; i < observers.size(); i++)
        delete observers[i];
}

void LogSimulation::open()
{
    std::cout << "otwieram plik do logowania" << std::endl;
    os = fopen(logfile, "w");
}

void LogSimulation::close()
{
    fclose(os);
}

void LogSimulation::readConfigFile()
{
    std::ifstream cfile;
    cfile.open(configfile);
    
    std::string line;
    while( std::getline (cfile, line) )
    {
        std::cout << line << '\n';
    }
    
    cfile.close();
}

void LogSimulation::registerObservers()
{
    readConfigFile();
    
    //SurfacePressure
            
    //std::string const s = "SurfacePressure";
            
    Observer* test_obj = ObserverFactory::createInstance("SurfacePressure","SP","%10.5f ");//("O1", "%10.5f ");
    
    //Observer* sp1 = new SurfacePressure("SurfacePressure1 ", "%10.5f ");
    //Observer* sp2 = new SurfacePressure("SurfacePressure2 ", "%6.3f ");
    //sp1->set_params(1,0.25);
    //sp2->set_params(1,0.123);
    
    //observers.push_back( sp1 );
    //observers.push_back( sp2 );
    std::cout << "jestem po" << std::endl;
}

void LogSimulation::printHeader()
{
    
    std::cout << "drukuje header" << std::endl;
    fprintf(os,"#");
    
    for(std::vector<Observer*>::iterator it = observers.begin(); it != observers.end(); ++it)
    {
        fprintf(os, (*it)->getName() );
    }
    
   fprintf(os,"\n");
}

void LogSimulation::dumpState(Box& box, std::vector<Cell>& cells, double rv, int simstep, int numV, int nbhandler)
{
    std::cout << "Dumpuje state. observers.size()="<< observers.size()<< std::endl;

    
    for(std::vector<Observer*>::iterator it = observers.begin(); it != observers.end(); ++it)
    {
       //std::cout << "jestem w petli" << std::endl;
       //std::cout << (*it)->observe(box, cells) << std::endl;
       fprintf(os, (*it)->getFormat(), (*it)->observe(box, cells) );
       
    }
    fprintf(os, "\n");
    
    //std::string FORMAT;
//    std::cout << "liczba observerow=" << observers.size() << std::endl;
//    for (int i = 0 ; i < observers.size(); i++)
//    {
//       fprintf(os, observers[i]->getFormat(), observers[i]->observe(box, cells) );
//    }
//
//    fprintf(os,"\n");
    
//    double pressure = SurfacePressure::calcPressure(box, cells, 0.0);
//    double volume_frac = VolumeFraction::calcVolumeFraction(box, cells, rv);
//    double area_coverage = WallCoverageFraction::wallsCoverage(box, cells, rv);
//    double mean_stress = AverageContactStress::caclContactStress(box, cells);
//    double tot_surface = TotalCellsArea::totalCellArea(cells);
//    double strain_energy = SurfaceStrainEnergy::calcSurfaceEnergy(cells);
//    double box_volume = box.getVolume(0.0);
//    double box_area = box.getArea(0.0);
//    int numofcells = cells.size();
//    double average_turgor = AverageTurgor::populationAverageTurgor(cells) ;
//    fprintf(os, "%i %i %8.6f %8.6f %8.6f %8.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", simstep, numofcells, box_volume, box_area, pressure, volume_frac, area_coverage, mean_stress, tot_surface, strain_energy, average_turgor);
//    fflush(os);
}