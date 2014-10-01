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
#include "utils/Logger.h"

#define STRCMP(a,b) (!strcmp(a,b))

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

        int getNumberOfCells()
        {
            return numberofCells;
        }
        void setBoxSize(double);
        int getTotalVertices();

    private:

        void (Simulator::*integrator)();
        void setIntegrator(void (Simulator::*functoall)());
        void diagnoseParams();
        void rebuildVerletLists();

        vector<Cell> cells;
        arguments params;
        int numberofCells;
        double dt;
        double a;
        double dp;
        double gamma;
        double R0;
        double Rc;
        double ttotal;
        double initcellmass;
        double verlet_r;
        int nsteps;
        int d;

        int logStep;
        int saveStep;
        int vlistStep;
        int boxStep;

        Box box;
        bool pbc;
        bool drawBox;

        ScriptBuilder sb;
        XyzTraj traj;

        utils::Logger simulator_logs;

};

#endif	/* SIMULATOR_H */

