#ifndef SIMULATOR_H
#define	SIMULATOR_H

#include <fstream>
#include <cstring>
#include <vector>
#include <stdio.h>      /* fprintf*/
#include <stdlib.h>     /* atoi,  strtod */
#include <assert.h>

#include "Box.h"
#include "Environment.h"
#include "Cell.h"
#include "arguments.h"
#include "exceptions/MaxSizeException.h"
#include "exceptions/DataException.h"
#include "exceptions/NotImplementedException.h"
#include "geometry/algorithms/SimpleTriangulation.h"
#include "geometry/algorithms/PlatonicTriangulatoin.h"
#include "utils/io/ScriptBuilder.h"
#include "utils/io/XyzTraj.h"
#include "utils/io/LogSimulation.h"
#include "utils/Logger.h"
#include "force/OsmoticForce.h"
#include "simulation/DomainList.h"
#include "simulation/Restarter.h"
#include "simulation/Energy.h"
#include "utils/nrutil.h"
//#include "cg/cg.h"

struct params_t
{
    int log_step;
    int save_step;
    int nsteps;
    int d;
    int nbhandler;
    int platotype;

    double E_cell;
    double nu;
    double th;
    double dt;
    double dp;
    double ddp;
    double volume_scale;
    double ttime;
    double r_vertex;
    bool draw_box;
    bool scale;
    bool dynamics;
    bool const_volume;
};

class Simulator
{
    public:
        Simulator(const arguments&);
        Simulator(const Simulator& orig) = delete;
        virtual ~Simulator();

        void simulate();
        void simulate(int);

        void initCells(int, double);
        void initCells(int, double, double);
        void initCells(int, double, double, char*, bool = false);

        void restart();

    private:

        void (Simulator::*integrator)();
        void shiftCell(const Vector3D&, int);
        void setIntegrator(void (Simulator::*functoall)());
        void setIntegrator(char*);
        void setTriangulator(char*);
        void diagnoseParams(arguments);
        void logParams();

        void pushCell(const Cell&);
        void addCell(double, char*);

        void calcForces();
        void integrate();
        void rebuildDomainsList();

        void integrateEuler();
        void heunMethod();
        void midpointRungeKutta();
        void gear_cp();

        int getTotalVertices();
        double getLengthScale();

        void updateCells();
        void update_neighbors_list();

        void set_min_force();
        bool check_min_force();
        bool check_const_volume();

        double MIN_FORCE_SQ = 0.0;
        double FORCE_FRAC = 0.01;

        int number_of_cells;
        std::string triangulator;
        params_t params;

        std::vector<Cell> cells;

        Box box;

        ScriptBuilder sb;
        XyzTraj traj;
        LogSimulation log_sim;

        DomainList domains;
        Restarter restarter;

        static utils::Logger simulator_logs;
        static ulong FORCE_EVALUATION_COUTER;

        // CONJUGATE GRADIENTS CODE
        void cg();
        void boost_cg();

        double func(double[]);
        void dfunc(double[], double[]);
        void frprmn(double p[], int n, double ftol, int* iter, double* fret);
        void linmin(double p[], double xi[], int n, double* fret );
        void dlinmin(double p[], double xi[], int n, double* fret );
        void mnbrak(double* ax, double* bx, double* cx, double* fa, double* fb, double* fc);
        double brent(double ax, double bx, double cx, double tol, double* xmins);
        double dbrent(double ax, double bx, double cx, double tol, double* xmin);
        double df1dim(double x);
        double f1dim(double x);
};

#endif	/* SIMULATOR_H */