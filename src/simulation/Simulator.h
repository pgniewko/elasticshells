#ifndef SIMULATOR_H
#define	SIMULATOR_H

#define MAX_CELLS 100
#define STRCMP(a,b) (!strcmp(a,b))

#include "../Cell.h"

using namespace std;
class Simulator {
public:
    Simulator();
    Simulator(const Simulator& orig);
    virtual ~Simulator();
    
    void integrateEuler();
    void integrateVv();
    void integrateDampedEuler();
    
    void integrate();
    void setIntegrator(char* token);
    
    void addCell(const Cell&);
    
    void calcForces();
    void move();
        
private:
    vector<Cell> cells;
    void (Simulator::*integrator)();
    int numberofCells;
    

};

#endif	/* SIMULATOR_H */

