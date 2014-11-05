#ifndef SIMULATOR_H
#define	SIMULATOR_H

#include <fstream>
#include <float.h>      /* DBL_MAX */
#include <cstring>
#include <vector>
#include <stdio.h>      /* fprintf*/
#include <stdlib.h>    /* atoi,  strtod */

#include "Box.h"
#include "Environment.h"
#include "Cell.h"
#include "arguments.h"
#include "exceptions/MaxSizeException.h"
#include "exceptions/DataException.h"
#include "exceptions/NotImplementedException.h"
#include "geometry/algorithms/SimpleTriangulation.h"
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
    int box_step;
    int vlist_step;
    int nsteps;
    int d;
    int nbhandler;
    
    double a;
    double dt;
    double dp;
    double visc;
    double k;
    double mass;
    double ttime;
    double r_cut;
    double r_bc;
    double verlet_r;
    bool draw_box;
};

class Simulator
{
    public:
        Simulator(const arguments&);
        Simulator(const Simulator& orig);
        virtual ~Simulator();

        void integrateEuler();
        void integrateVv();
        void heunMethod();
        void midpointRungeKutta();


        void simulate();
        void simulate(int);
        void integrate();
        void setIntegrator(char* token);

        void addCell(const Cell&);
        void addCell();
        void addCell(double);

        void initCells(int, double);
        void initCells(int, double, double);

        void calcForces();

        void moveCell(const Vector3D&, int);
        void addCellVel(const Vector3D&, int);

        int getNumberOfCells();
        void setBoxSize(double);
        int getTotalVertices();

    private:

        void (Simulator::*integrator)();
        void setIntegrator(void (Simulator::*functoall)());
        void diagnoseParams(arguments);
        void logParams();
        void rebuildVerletLists();
        void rebuildDomainsList();
        double getMaxScale();

        int numberofCells;
        params_t params;
        
        std::vector<Cell> cells;
        Box box;

        ScriptBuilder sb;
        XyzTraj traj;
        LogSimulation logsim;
        utils::Logger simulator_logs;
        
        DomainList domains;
};

#endif	/* SIMULATOR_H */

