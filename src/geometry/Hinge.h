#ifndef HINGE_H
#define	HINGE_H

#include <iostream>
#include <iomanip>

#include "Environment.h"
#include "geometry/Triangle.h"
#include "geometry/Vertex.h"
#include "geometry/Vector3D.h"

class Hinge
{
        friend class Restarter;
        friend class Simulator;
        friend class Observer;
    public:
        Hinge();
        explicit Hinge(int, int, int, int);
        Hinge(const Hinge& orig);
        virtual ~Hinge();

        void set_d(const double&, const double&, const double&);
        double calc_radius_of_curvature(std::vector<Vertex>&) const;
        void set_theta_zero(const std::vector<Vertex>&);

        void set_id(int);
        int get_id() const;

        friend bool operator== (Hinge&, Hinge&);
        friend std::ostream& operator<< (std::ostream&, const Hinge&);

    private:
        double D = 0.0;
        double sin_theta_0 = 0.0;
        double theta_0 = 0.0;
        int x1 = -1;
        int x2 = -1;
        int x3 = -1;
        int x4 = -1;
        int my_id = 1;
        double calc_theta(const std::vector<Vertex>&) const;
        double calc_sin_theta(const std::vector<Vertex>&) const;
};

inline bool operator== (Hinge& bs1, Hinge& bs2)
{
    return (bs1.x1 == bs2.x1 && bs1.x2 == bs2.x2 && bs1.x3 == bs2.x3 && bs1.x4 == bs2.x4);
}

#endif	/* HINGE_H */

