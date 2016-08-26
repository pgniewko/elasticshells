#ifndef BENDINGHINGE_H
#define	BENDINGHINGE_H

#include <iostream>
#include <iomanip>

#include "Environment.h"
#include "geometry/Triangle.h"
#include "geometry/Vertex.h"
#include "geometry/Vector3D.h"

class BendingHinge
{
        friend class Restarter;
    public:
        BendingHinge();
        explicit BendingHinge(int, int, int, int);
        BendingHinge(const BendingHinge& orig);
        virtual ~BendingHinge();

        void setD(const double&, const double&, const double&);
        void calcBendingForces(Vertex[]) const;
        double calcRadiusOfCurvature(Vertex[]) const;
        void setThetaZero(const Vertex[]);
        double calcBendingEnergy(const Vertex[]) const;

        void setId(int);
        int getId() const;

        friend bool operator== (BendingHinge&, BendingHinge&);
        friend std::ostream& operator<< (std::ostream&, const BendingHinge&);

    private:
        double D = 0.0;
        double sinTheta0 = 0.0;
        double theta0 = 0.0;
        int x1 = -1;
        int x2 = -1;
        int x3 = -1;
        int x4 = -1;
        int myid = 1;
        double calcTheta(const Vertex[]) const;
        double calcSinTheta(const Vertex[]) const;
};

inline bool operator== (BendingHinge& bs1, BendingHinge& bs2)
{
    return (bs1.x1 == bs2.x1 && bs1.x2 == bs2.x2 && bs1.x3 == bs2.x3 && bs1.x4 == bs2.x4);
}

#endif	/* BENDINGHINGE_H */
