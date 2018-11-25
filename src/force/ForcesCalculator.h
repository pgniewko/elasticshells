#ifndef FORCESCALCULATOR_H
#define FORCESCALCULATOR_H

#include "Environment.h"
#include "geometry/Vector3D.h"
#include "geometry/Tetrahedron.h"
#include <domain_list_t.h>
#include "simulation/elements.h"


class ForcesCalculator {
public:
    ForcesCalculator(int, bool, bool);
    ForcesCalculator(const ForcesCalculator& orig);
    virtual ~ForcesCalculator();
    
    void calculate_forces(const std::vector<double>& xyz, 
                          std::vector<double>& forces,
                          const std::vector<element>& elements,
                          const std::vector<hinge>& hinges,
                          const std::vector<object_map> vs_map,
                          const std::vector<double> turgors,
                          const int num_shells) const;
    
private:
    int m;
    bool pbc;
    bool bending;
    domain_list_t dl;
    
    void evaluate_elements(const std::vector<double>& xyz,
                           std::vector<double>& forces,
                           const std::vector<element>& elements) const;
    
    void evaluate_hinges( const std::vector<double>& xyz,
                          std::vector<double>& forces,
                          const std::vector<hinge>& hinges) const;
    
    void evaluate_pressure(const std::vector<double>& xyz, 
                           std::vector<double>& forces,
                           const std::vector<element>& elements,
                           const std::vector<object_map> vs_map,
                           const std::vector<double> turgors,
                           const int num_shells) const;
    
    double calculate_theta(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&) const;
    Vector3D calculate_dV(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&) const;
};

#endif /* FORCESCALCULATOR_H */

