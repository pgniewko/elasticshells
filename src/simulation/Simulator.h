#ifndef SIMULATOR_H
#define	SIMULATOR_H

#include <fstream>
#include <cstring>
#include <vector>
#include <stdio.h>      /* fprintf*/
#include <stdlib.h>     /* atoi,  strtod */
#include <assert.h>
#include <set>
#include <unordered_map>

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
#include "utils/io/XyzTraj.h"
#include "utils/io/LogSimulation.h"
#include "utils/Logger.h"
#include "force/OsmoticForce.h"
#include "simulation/Restarter.h"
#include "simulation/Energy.h"
#include "simulation/Packer.h"
#include "utils/nrutil.h"
#include "../Timer.h"
#include "elements.h"

#include "integrators/Integrator.h"
#include "../force/ForcesCalculator.h"

#define MAX_M 200

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

class Integrator;

class Simulator
{
        friend class Integrator;

    public:
        Simulator(const arguments&);
        Simulator(const Simulator& orig) = delete;
        virtual ~Simulator();

        void simulate(int);
        void restart();
        void analyze();
        
        void init_shells(int, double, double, bool = false);

    private:

        void shiftShell(const Vector3D&, int);
        void setTriangulator(char*);
        void diagnoseParams(arguments);
        void logParams();

        void pushShell(const Shell&);
        void addShell(double);

        void calcForces();
        void integrate();

        int getTotalVertices();
        double getLengthScale(double = 0.0);

        void recalculate_mass_centers();
        void set_min_force();
        bool check_min_force();
        bool check_const_volume();

        void saveTurgors();

        void recenterShells();

        double volumeFraction();
        
        void create_shells_image();
        
        void copy_shells_data();
        void copy_back_shells_data();
        
        int estimate_m();

        Integrator* integrator;

        int number_of_shells;
        std::string triangulator;
        params_t params;

        std::vector<Shell> shells;
        std::vector<double> xyz;
        std::vector<double> forces;
        
        int vertex_number = 0;
        int hinge_number = 0 ;
        int element_number = 0;
        
        std::vector<element> elements;
        std::vector<hinge> hinges;
        std::vector<object_map> vs_map;
        std::map<object_map, int> inv_vs_map;
        
        std::vector<object_map> ts_map;
        std::map<object_map, int> inv_ts_map;
        
        std::vector<object_map> hs_map;
        std::map<object_map, int> inv_hs_map;
        
        std::vector<double> turgors;
        
        std::vector<std::vector<int> > graph;

        Box box;

        XyzTraj traj;
        LogSimulation log_sim;

        Restarter restarter;
        
        ForcesCalculator fc;
        double MIN_FORCE = 0.0;
        double FORCE_FRAC = 0.0001;


        static utils::Logger simulator_logs;
        static unsigned long FORCE_EVALUATION_COUTER;

        static bool RESTART_FLAG;
};

#endif	/* SIMULATOR_H */
