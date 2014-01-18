#include "RepulsiveForce.h"

RepulsiveForce::RepulsiveForce(double a, double r_cut) {
    this->a = a;
    this->r_cut = r_cut;
    cout << a << " " << r_cut << endl;
}

RepulsiveForce::RepulsiveForce(const RepulsiveForce& orig) {
}

RepulsiveForce::~RepulsiveForce() {
}

inline double cforce(double r, double rcut)
{
    return (1.0 -  r / rcut);
}

Vector3D RepulsiveForce::eval(Vector3D r_kl, Cell cell_k, Cell cell_l)
{
    //return Vector3D();
    //cout << a << r_cut << endl;
    double R_kl = r_kl.length();
    //k
    //Vector3D U_k = cell_k.p / cell_k.mass ;
    //l
    //Vector3D U_l = cell_l.p / cell_l.mass;
    //geometry
    Vector3D e_kl(r_kl / R_kl);




    Vector3D Fc_kl = a * cforce(R_kl, r_cut) * e_kl;
    return Fc_kl;
}

inline double RepulsiveForce::magn(Vector3D r_kl, Cell cell_k, Cell cell_l)
{
    //return Vector3D();
    //cout << a << r_cut << endl;
    double R_kl = r_kl.length();
    //k
    //Vector3D U_k = cell_k.p / cell_k.mass ;
    //l
    //Vector3D U_l = cell_l.p / cell_l.mass;
    //geometry
    //Vector3D e_kl(r_kl / R_kl);




    double Fc_kl = a * cforce(R_kl, r_cut);// * e_kl;
    return Fc_kl;
}