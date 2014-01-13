#ifndef SIMULATOR_H
#define	SIMULATOR_H

#include <vector>
#include <fstream>
#include <cmath>
#include "Vector3D.h"
#include "Pairlist.h"
#include "random.h"
#include "constants.h"
#include "arguments.h"

using namespace std;

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
        //Simulator(const Parameters& par, const Distance& _domain);
        Simulator(const arguments& par, const Distance& _domain);
//        virtual ~Simulator();

        void integrate_euler();
        void integrate_DPD_VV();
        void integrate_trotter();
        void state(double& T, Vector3D& P) const;
        int read_pos(const char* filename);
        void write_pos(const char* filename, bool wrap = 1);
        arguments params;

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

#endif	/* SIMULATOR_H */

