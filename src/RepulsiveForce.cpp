#include "RepulsiveForce.h"

RepulsiveForce::RepulsiveForce(double a, double r_cut) {
    this->a = a;
    this->r_cut = r_cut;
}

//RepulsiveForce::RepulsiveForce(const RepulsiveForce& orig) {
//}

//RepulsiveForce::~RepulsiveForce() {
//}

inline double cforce(double r, double rcut)
{
    return (1.0 -  r / rcut);
}

Vector3D RepulsiveForce::eval(Vector3D r_kl, Cell cell_k, Cell cell_l)
{
    double R_kl = r_kl.length();
    Vector3D e_kl(r_kl / R_kl);
    return a * cforce(R_kl, r_cut) * e_kl;
}

inline double RepulsiveForce::magn(Vector3D r_kl, Cell cell_k, Cell cell_l)
{
    double R_kl = r_kl.length();
    return a * cforce(R_kl, r_cut);
}