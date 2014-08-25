#ifndef SIMULATOR_H
#define	SIMULATOR_H

#define MAX_CELLS 1000

#include <cstring>
#include <vector>

#include "../Cell.h"
#include "../arguments.h"
#include "../geometry/algorithms/SimpleTriangulation.h"

using namespace std;

#define STRCMP(a,b) (!strcmp(a,b))

class Simulator {
public:
    Simulator(const arguments&);
    Simulator(const Simulator& orig);
    virtual ~Simulator();
    
    void integrateEuler();
    void integrateVv();
    void integrateDampedEuler();
    
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
        
private:
    
    int ncells;
    
    vector<Cell> cells;
    void (Simulator::*integrator)();
    void setIntegrator(void (Simulator::*functoall)());
    int numberofCells;
    
    double dt;
    double a;
    double dp;
    double gamma;
    double R0;
    double Rc;
    double ttotal;
    int nsteps;
    int d;
    

};

#endif	/* SIMULATOR_H */

