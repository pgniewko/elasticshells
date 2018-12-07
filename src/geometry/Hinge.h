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

        void setD(const double&, const double&, const double&);
        double calcRadiusOfCurvature(std::vector<Vertex>&) const;
        void setThetaZero(const std::vector<Vertex>&);

        void setId(int);
        int getId() const;

        friend bool operator== (Hinge&, Hinge&);
        friend std::ostream& operator<< (std::ostream&, const Hinge&);

    private:
        double D = 0.0;
        double sinTheta0 = 0.0;
        double theta0 = 0.0;
        int x1 = -1;
        int x2 = -1;
        int x3 = -1;
        int x4 = -1;
        int myid = 1;
        double calcTheta(const std::vector<Vertex>&) const;
        double calcSinTheta(const std::vector<Vertex>&) const;
};

inline bool operator== (Hinge& bs1, Hinge& bs2)
{
    return (bs1.x1 == bs2.x1 && bs1.x2 == bs2.x2 && bs1.x3 == bs2.x3 && bs1.x4 == bs2.x4);
}

#endif	/* HINGE_H */

