#ifndef SIMULATOR_H
#define	SIMULATOR_H

#include <vector>
#include <fstream>
#include <cmath>
#include <cstring>

#include "constants.h"
#include "arguments.h"
#include "random.h"
#include "Cell.h"
#include "Force.h"
#include "Pairlist.h"
#include "RepulsiveForce.h"
#include "Vector3D.h"

using namespace std;

class Simulator
{
    public:
        Simulator(const arguments& par, const Box& _domain);
        virtual ~Simulator();

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
        void set_integrator(char* token);
        void set_forces();
        arguments params;


    private:
        int np;

        Box box;

        Pairlist pl;

        Force* cForce;
        void (Simulator::*integrator)();

        vector<Pair> pairs;
        vector<Cell> cells;
        vector<Vector3D> forces;

        double compute_forces(double dt);
        double compute_forces_trotter(double dt, int order = 0);
        void set_integrator(void (Simulator::*functoall)());
};

#endif	/* SIMULATOR_H */