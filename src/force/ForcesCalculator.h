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
                          const int num_shells, 
                          const double rv, const double E, const double nu,
                          const double Eb, const double nub);
    
private:
    int m;
    bool pbc;
    bool bending;
    domain_list_t dl;
    
    void evaluate_elements(const std::vector<double>& xyz,
                           std::vector<double>& forces,
                           const std::vector<element>& elements) const;
    
    void evaluate_hinges(const std::vector<double>& xyz,
                         std::vector<double>& forces,
                         const std::vector<hinge>& hinges) const;
    
    void evaluate_pressure(const std::vector<double>& xyz, 
                           std::vector<double>& forces,
                           const std::vector<element>& elements,
                           const std::vector<object_map> vs_map,
                           const std::vector<double> turgors,
                           const int num_shells) const;
    
    void evaluate_nonbonded(const std::vector<double>& xyz, 
                            std::vector<double>& forces,
                            const std::vector<element>& elements,
                            const double rv, const double E, const double nu);
    
    void evaluate_box(const std::vector<double>& xyz, 
                            std::vector<double>& forces,
                            const double rv, const double E, const double nu, 
                            const double Eb, const double nub);
    
    double calculate_theta(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&) const;
    Vector3D calculate_dV(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&) const;
    void distance(Vector3D& dkj, const Vector3D& vj, const Vector3D& vk) const;
};

#endif /* FORCESCALCULATOR_H */

