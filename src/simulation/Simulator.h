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
#include "Shell.h"
#include "arguments.h"
#include "exceptions/MaxSizeException.h"
#include "exceptions/DataException.h"
#include "exceptions/NotImplementedException.h"
#include "exceptions/NotAllowedException.h"
#include "geometry/algorithms/SimpleTriangulation.h"
#include "geometry/algorithms/PlatonicTriangulatoin.h"
#include "geometry/algorithms/RandomTriangulation.h"
#include "utils/io/ScriptBuilder.h"
#include "utils/io/XyzTraj.h"
#include "utils/io/LogSimulation.h"
#include "utils/Logger.h"
#include "force/OsmoticForce.h"
#include "simulation/DomainList.h"
#include "simulation/Restarter.h"
#include "simulation/Energy.h"
#include "simulation/Packer.h"
#include "utils/nrutil.h"
#include "../Timer.h"

#include "integrators/Integrator.h"

struct params_t
{
    int log_step;
    int nsteps;
    int d;
    int nbhandler;
    int platotype;

    double E_shell;
    double nu;
    double th;
    double dt;
    double dp;
    double ddp;
    double ttime;
    double r_vertex;
    bool draw_box;
    bool const_volume;
};

struct element
{
    int ia, ib, ic;
    double an[3];
    double L2[3];
    double ki[3];
    double ci[3];
};

struct hinge
{
    int x1, x2, x3, x4;
    double D;
    double sinTheta0;
    double theta0;
};

class Integrator;

class Simulator
{
        friend class Integrator;

    public:
        Simulator(const arguments&);
        Simulator(const Simulator& orig) = delete;
        virtual ~Simulator();

        void simulate();
        void simulate(int);

        void initShells(int, double, double, bool = false);

        void restart();
        void analyze();

    private:

        void shiftShell(const Vector3D&, int);
        void setTriangulator(char*);
        void diagnoseParams(arguments);
        void logParams();

        void pushShell(const Shell&);
        void addShell(double);

        void calcForces();
        void integrate();
        void rebuildDomainsList();


        int getTotalVertices();
        double getLengthScale(double = 0.0);

        void updateShells();
        void update_neighbors_list();

        void set_min_force();
        bool check_min_force();
        bool check_const_volume();

        void saveTurgors();

        void recenterShells();

        double volumeFraction();

        Integrator* integrator;

        int number_of_shells;
        std::string triangulator;
        params_t params;

        std::vector<Shell> shells;
        
        // ****
        std::vector<double> xyz;
        std::vector<double> forces;
        
        int vertex_number = 0;
        int hinge_number = 0 ;
        int element_number = 0;
        
        std::vector<element> elements;
        std::vector<hinge> hinges;
        /// ************
        

        Box box;

        ScriptBuilder sb;
        XyzTraj traj;
        LogSimulation log_sim;

        DomainList domains;
        Restarter restarter;

        double MIN_FORCE_SQ = 0.0;
        double FORCE_FRAC = 0.0001;


        static utils::Logger simulator_logs;
        static unsigned long FORCE_EVALUATION_COUTER;

        static bool RESTART_FLAG;

};

#endif	/* SIMULATOR_H */
