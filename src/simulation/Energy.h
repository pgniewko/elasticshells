#ifndef ENERGY_H
#define	ENERGY_H

#include "Environment.h"
#include "simulation/Box.h"
#include "Cell.h"

class Energy {
public:
    Energy();
    Energy(const Energy& orig);
    virtual ~Energy();
    
    static double calcTotalEnergy(const std::vector<Cell>&, const Box&);
    static double calcMembraneEnergy(const std::vector<Cell>&, const Box&);
    static double calcOsmoticEnergy(const std::vector<Cell>&);
    static double calcContactEnergy(const std::vector<Cell>&, const Box&);
    
    static double calcBendingEnergy(const std::vector<Cell>&);
    static double calcStretchEnergy(const std::vector<Cell>&, const Box&, char*);
    static double calcCellCellEnergy(const std::vector<Cell>&, const Box&);
    static double calcCellBoxEnergy(const std::vector<Cell>&, const Box&);
    
    
    static ulong ENERGY_EVALUATION_COUNTER;
    
private:

    static utils::Logger energy_logs;
};

#endif	/* ENERGY_H */

