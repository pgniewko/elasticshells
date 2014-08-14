#ifndef SIMULATOR_H
#define	SIMULATOR_H

#define MAX_CELLS 100

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
        
private:
    void (Simulator::*integrator)();
    

};

#endif	/* SIMULATOR_H */

