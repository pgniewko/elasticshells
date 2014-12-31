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

void LogSimulation::dumpState(Box& box, std::vector<Cell>& cells, double rv, int simstep, int numV, int nbhandler)
{
    double pressure = SurfacePressure::calcPressure(box, cells, rv);
    double volumeFrac = VolumeFraction::caclVolumeFraction(box, cells, rv);
    double cellsVolume = VolumeFraction::caclCellsVolume(cells);
    double boxVolume = box.getVolume(rv);
    int numofcells = cells.size();
    fprintf(os, "%i %i %i %8.2f %8.2f %10.8f %8.6f\n", simstep, numofcells, numV, boxVolume, cellsVolume, volumeFrac, pressure);
    //double q6 = QL::calcQl(cells[0], 4, 5.0);
    //fprintf(os, "%i %i %i %8.2f %8.2f %10.8f %8.6f %8.6f\n", simstep, numofcells, numV, boxVolume, cellsVolume, volumeFrac, pressure, q6);
}