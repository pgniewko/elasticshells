#include "LogSimulation.h"

LogSimulation::LogSimulation(char* lf) 
{
    logfile = lf;
}

LogSimulation::LogSimulation(const LogSimulation& orig) : logfile(orig.logfile) {}

LogSimulation::~LogSimulation() {}

void LogSimulation::open()
{
    os = fopen(logfile, "w");
}

void LogSimulation::close()
{
    fclose(os);
}

void LogSimulation::dumpState(Box& box, std::vector<Cell>& cells, int simstep, int numV, int nbhandler)
{
    
    //double pressure = SurfacePressure::calcPressure(box, cells);
    //double volumeFrac = VolumeFraction::caclVolumeFraction(box, cells);
    //double cellsVolume = VolumeFraction::caclCellsVolume(cells);
    //double volPressure = VolumePressure::calcPressure(box, cells);
    //double boxVolume = box.getVolume();
    //int numofcells = cells.size();
    //fprintf(os, "%i %i %i %8.2f %8.2f %10.8f %8.6f %8.6f\n", simstep, numofcells, numV, boxVolume, cellsVolume, volumeFrac, pressure, volPressure);
}