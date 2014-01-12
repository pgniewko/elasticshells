#ifndef SIMULATOR_H
#define	SIMULATOR_H

#include <vector>
#include <fstream>
#include <cmath>
#include "Vector3D.h"
#include "Pairlist.h"
#include "random.h"
#include "constants.h"

using namespace std;

struct Parameters
{
    int Nparticles;
    double L;
    int Niter;
    double rcut;
    double pairdist;
    double tau;
    char* infname;
    double a;
    double gamma;
    double sigma;
    double mass;
};


class Fields {
    public:
        Vector3D P;
        Fields() : P(0.0, 0.0, 0.0) {}
};

class Pair {
    public:
        int k, l;
        Pair(int _k, int _l) : k(_k), l(_l) {}
};

class Simulator {
    public:
//    Simulator(Vector3D* pp, Fields* fp, int _np, const Distance& _domain, const Parameters& par);
        Simulator(const Parameters& par, const Distance& _domain);
        virtual ~Simulator();

        void integrate_euler();
        void integrate_DPD_VV();
        void integrate_trotter();
        void state(double& T, Vector3D& P) const;
        int read_pos(const char* filename);
        void write_pos(const char* filename, bool wrap = 1);
        Parameters Par;

    private:
        Distance domain;
        Pairlist pl;
        vector<Pair> pairs;
        vector<Fields> fp;
        vector<Vector3D> r;
        int np;
        double compute_forces(double dt, vector<Fields>& dfp, int NC = 1);
        double compute_forces_trotter(double dt, int order = 0);
};


//class Simulator {
//public:
//    Simulator();
//    Simulator(const Simulator& orig);
//    virtual ~Simulator();

//    int generate_pos();
//    int read_pos(const char* filename);
//    void write_pos(const char* filename, bool wrap=1);
//private:

//};

#endif	/* SIMULATOR_H */

