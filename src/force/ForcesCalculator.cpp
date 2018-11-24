#include "ForcesCalculator.h"

ForcesCalculator::ForcesCalculator(int m_, bool pbc_, bool bend) : m(m_), pbc(pbc_), bending(bend), dl(m_,pbc_) 
{
    
}

ForcesCalculator::ForcesCalculator(const ForcesCalculator& orig) : m(orig.m), pbc(orig.pbc), bending(orig.bending), dl(orig.m,orig.pbc) {}

ForcesCalculator::~ForcesCalculator() {
}

void ForcesCalculator::calculate_forces(std::vector<double>& xyz, std::vector<double>& forces,
        std::vector<element>& elements, std::vector<hinge>& hinges, std::vector<object_map> vs_map)
{
    
    
    // ITERATE OVER ELEMENTS
    
    // ITERATE OVER HINGES
    
    // CALCULATE MASS CENTERS
    // PRESSURE FORCES
    
    // END WITH NON-BONDED
    
    
}

