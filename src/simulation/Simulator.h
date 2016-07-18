#ifndef SIMULATOR_H
#define	SIMULATOR_H

#include <fstream>
#include <cstring>
#include <vector>
#include <stdio.h>      /* fprintf*/
#include <stdlib.h>     /* atoi,  strtod */

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
        void initCells(int, double, double, char*);

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

        static utils::Logger simulator_logs;
};

#endif	/* SIMULATOR_H */
