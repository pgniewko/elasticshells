#ifndef FORCESCALCULATOR_H
#define FORCESCALCULATOR_H

#include "Environment.h"
#include "geometry/Vector3D.h"
#include "geometry/Tetrahedron.h"
#include <domain_list_t.h>
#include "simulation/elements.h"
#include "Shell.h"
#include <set>
#include <unordered_map>


class ForcesCalculator
{
    public:
        ForcesCalculator();
        ForcesCalculator(int, bool, bool);
        ForcesCalculator(const ForcesCalculator& orig);
        virtual ~ForcesCalculator();

        void calculate_forces(const std::vector<double>& xyz,
                              std::vector<double>& forces,
                              const std::vector<element>& elements,
                              const std::vector<hinge>& hinges,
                              const std::vector<object_map>& vs_map,
                              const std::vector<std::vector<int> >& graph_,
                              const std::vector<double>& turgors,
                              const int num_shells,
                              const double rv, const double E, const double nu,
                              const double Eb, const double nub);
        void set_dl_dims(const double, const double, const int);

        void reset_dl(int, const Box&);

        pairs_t contacts_list(const std::vector<double>& xyz,
                              const std::vector<std::vector<int> >& graph,
                              const std::vector<object_map>& vs_map,
                              const int num_shells,
                              const double rv);

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
                               const std::vector<object_map>& vs_map,
                               const std::vector<double>& turgors,
                               const int num_shells) const;

        void evaluate_nonbonded(const std::vector<double>& xyz,
                                std::vector<double>& forces,
                                const std::vector<std::vector<int> >& graph,
                                const double rv, const double E, const double nu);

        void evaluate_box(const std::vector<double>& xyz,
                          std::vector<double>& forces,
                          const double rv, const double E, const double nu,
                          const double Eb, const double nub);

        void zero_forces(std::vector<double>& forces) const;

        double calculate_theta(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&) const;
        Vector3D calculate_dV(const Vector3D&, const Vector3D&, const Vector3D&, const Vector3D&) const;
        void distance(Vector3D& dkj, const Vector3D& vj, const Vector3D& vk) const;

        bool is_bonded(int, int, const std::vector<std::vector<int> >&);
};

#endif /* FORCESCALCULATOR_H */

