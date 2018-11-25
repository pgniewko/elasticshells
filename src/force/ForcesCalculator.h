#ifndef FORCESCALCULATOR_H
#define FORCESCALCULATOR_H

//#include <domain_list_t.h>
//#include "simulation/elements.h"


class ForcesCalculator {
public:
    ForcesCalculator(int, bool, bool);
    ForcesCalculator(const ForcesCalculator& orig);
    virtual ~ForcesCalculator();
    
//    void calculate_forces(std::vector<double>&,std::vector<double>&,std::vector<element>&,std::vector<hinge>&,std::vector<object_map>);
    
private:
    int m;
    bool pbc;
    bool bending;
//    domain_list_t dl;

};

#endif /* FORCESCALCULATOR_H */

