#ifndef REPULSIVEFORCE_H
#define	REPULSIVEFORCE_H

#include "Cell.h"
#include "Force.h"

using namespace std;

class RepulsiveForce : public Force
{
    public:
        RepulsiveForce() : a(25.0), r_cut(1.0) {};
        RepulsiveForce(double a, double r_cut);
        RepulsiveForce(const RepulsiveForce& orig) : a(orig.a), r_cut(orig.r_cut) {};
//    ~RepulsiveForce();
        Vector3D eval(Vector3D r_kl, Cell cell_k, Cell cell_l);
        inline double magn(Vector3D r_kl, Cell cell_k, Cell cell_l);
    private:
        double a;
        double r_cut;

};

#endif	/* REPULSIVEFORCE_H */
