#ifndef SIMULATOR_H
#define	SIMULATOR_H

//#define MAX_CELLS 100

#include <fstream>
#include <float.h>      /* DBL_MAX */
#include <cstring>
#include <vector>

#include "Box.h"
#include "Environment.h"
#include "Cell.h"
#include "arguments.h"
#include "exceptions/MaxSizeException.h"
#include "exceptions/DataException.h"
#include "exceptions/NotImplementedException.h"
#include "geometry/algorithms/SimpleTriangulation.h"


using namespace std;

#define STRCMP(a,b) (!strcmp(a,b))

class Simulator {
public:
    Simulator(const arguments&);
    Simulator(const Simulator& orig);
    virtual ~Simulator();
    
    void integrateEuler();
    void integrateVv();
    void heunMethod();
    void midpointRungeKutta();

    
    void simulate();
    void simulate(int);
    void integrate();
    void setIntegrator(char* token);
    
    void addCell(const Cell&);
    void addCell();
    
    void calcForces();
    void move();
    
    arguments params;
    
    void moveCell(const Vector3D&, int);
    void addCellVel(const Vector3D&, int);
    
    int getNumberOfCells() {return numberofCells;}
    
    void saveCellsState();
    void saveCellsState(const char*);
    void renderScript(bool box=false);
    void saveSurfaceScript();
    
    void setBoxSize(double);
    
    void printCell(int);
    
    int getTotalVertices();
    
    void printBox(ofstream&);
    
    Cell getCell(int);
        
private:
    
    vector<Cell> cells;
    void (Simulator::*integrator)();
    void setIntegrator(void (Simulator::*functoall)());
    void diagnoseParams();
    void rebuildVerletLists();
    int numberofCells;
    
    double dt;
    double a;
    double dp;
    double gamma;
    double R0;
    double Rc;
    double ttotal;
    double initcellmass;
    double verlet_r;
    int nsteps;
    int d;
    //double bs;
    
    int logStep;
    int saveStep;
    int vlistStep;
    int boxStep;
    
    Box box;
    //int boxSize;
    
    char* trajfile;
    char* script;
    char* surfaceScript;
    
    bool pbc;
    bool drawBox;
};

#endif	/* SIMULATOR_H */

