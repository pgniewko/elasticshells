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

void LogSimulation::registerObservers()
{

}

void LogSimulation::dumpState(Box& box, std::vector<Cell>& cells, double rv, int simstep, int numV, int nbhandler)
{
    double pressure = SurfacePressure::calcPressure(box, cells, rv);
    double volume_frac = VolumeFraction::caclVolumeFraction(box, cells, rv);
    double area_coverage = WallCoverageFraction::wallsCoverage(box, cells, rv);
    double mean_stress = AverageContactStress::caclContactStress(box, cells);
    double tot_surface = TotalCellsArea::totalCellArea(cells);
    double strain_energy = SurfaceStrainEnergy::calcSurfaceEnergy(cells);
    double box_volume = box.getVolume(0.0);
    double box_area = box.getArea(0.0);
    int numofcells = cells.size();
    double average_turgor = AverageTurgor::populationAverageTurgor(cells) ;
    fprintf(os, "%i %i %8.6f %8.6f %8.6f  %8.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", simstep, numofcells, box_volume, box_area, pressure, volume_frac, area_coverage, mean_stress, tot_surface, strain_energy, average_turgor);
    fflush(os);
    //double forces = SurfaceForce::calcForces(box, cells, rv);
    //double cellsVolume = VolumeFraction::caclCellsVolume(cells);
    //double boxVolume = box.getVolume(rv);
    //int numofcells = cells.size();
    //fprintf(os, "%i %i %i %8.2f %8.2f %10.8f %8.6f\n", simstep, numofcells, numV, boxVolume, cellsVolume, volumeFrac, pressure);
    //double ex = box.getXEdge(rv);
    //double ey = box.getYEdge(rv);
    //double ez = box.getZEdge(rv);
    //double q6 = QL::calcQl( cells, 6, box.getX() );
    //fprintf(os, "%9i %3i %5i %8.2f %8.2f %10.8f %8.6f %8.6f %8.6f %5.3f %5.3f %5.3f\n", simstep, numofcells, numV, boxVolume, cellsVolume, volumeFrac, pressure, q6, forces, ex, ey, ez);
}