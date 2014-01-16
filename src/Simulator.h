#ifndef SIMULATOR_H
#define	SIMULATOR_H

#include <vector>
#include <fstream>
#include <cmath>
#include <cstring>
#include "Vector3D.h"
#include "Pairlist.h"
#include "random.h"
#include "constants.h"
#include "arguments.h"
#include "Cell.h"
#include "YeastCell.h"

#include <error.h>     /* error */

using namespace std;

//class Fields {
//    public:
//        Vector3D P;
//        Fields() : P(0.0, 0.0, 0.0) {}
//};

class Pair
{
    public:
        int k, l;
        Pair(int _k, int _l) : k(_k), l(_l) {}
};


class Simulator
{
    public:
//    Simulator(Vector3D* pp, Fields* fp, int _np, const Distance& _domain, const Parameters& par);
        //Simulator(const Parameters& par, const Distance& _domain);
        Simulator(const arguments& par, const Box& _domain);
//        virtual ~Simulator();

        void integrate_euler();
        void integrate_DPD_VV();
        void integrate_trotter();
        void state(double& T, Vector3D& P) const;
        int initialize_pos(const char* filename, int np);
        int read_pos(const char* filename);
        int random_pos(int np);
        void write_pos(const char* filename, bool wrap = 1);
        void write_pos_traj(const char* filename, bool wrap = 1);
        void integrate();
        void set_integrator(void (Simulator::*functoall)() );
        void set_integrator(char* token);
        arguments params;

    private:
        Box box;
        Pairlist pl;
        vector<Pair> pairs;
        //vector<Fields> fp;
        //vector<Vector3D> r;

        vector<Cell> cells;

        int np;

        void (Simulator::*integrator)();
        double compute_forces(double dt, vector<Vector3D>& dfp, int NC = 1);
        double compute_forces_trotter(double dt, int order = 0);
};

#endif	/* SIMULATOR_H */

