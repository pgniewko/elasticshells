/*
 * Author : Pawel Gniewek (UC Berkeley)
 * Email  : pawel.gniewek@berkeley.edu
 * License: BSD
 */

#include <iostream>    /* cout, cin */
#include <fstream>     /* ios, ofstream*/
#include <stdio.h>     /* printf, fgets */
#include <argp.h>      /* argp_parse */
#include <stdlib.h>    /* atoi,  strtod */ 
#include <math.h>      /* log, sqrt */
#include <list>


#include "src/geometry/algorithms/SimpleTriangulation.h"
#include "src/geometry/Vector3D.h"
#include "src/geometry/Triangle.h"
#include "src/geometry/Vertex.h"
#include "src/arguments.h"
#include "src/random.h"
#include "src/Timer.h"
#include "src/Cell.h"
#include "src/force/HookeanForce.h"

using namespace std;

const char* argp_program_version = "biofilm 0.1.0";
const char* argp_program_bug_address = "<pawel.gniewek@berkeley.edu>";

static char doc[] =
    "General information about the software\
    \vThis part of the documentation comes *after* the options;";

static char args_doc[] = "\t\t\t(1st form)\n[STRING...]\t\t(2nd form)";

#define OPT_ABORT  1            /* --abort ^A*/

static struct argp_option options[] =
{
    {0,            0, 0, 0, "I/O Options:", 1},
    {"verbose",  'v', "INT", OPTION_ARG_OPTIONAL, "Produce verbose output for values [default: 1]" },
    {"debug",    'd', 0, OPTION_ALIAS},
    {"quiet",    'q', "INT", OPTION_ARG_OPTIONAL, "Don't produce any output [default: 1]" },
    {"silent",   's', 0, OPTION_ALIAS},
    {"input",    'i', "FILE",  0, "Input from FILE [default: ...]" },
    {"output",   'o', "FILE",  0, "Output to FILE instead of standard output [default: ... ]" },
    {"log",      'l', "FILE",  0, "Print log to FILE instead of standard output [default: ... ]" },
    {"xyz",      't', "FILE",  0, "Print trajectory to FILE [default: ... ]" },
    {"abort", OPT_ABORT, 0, 0, "Abort before showing any output"},

    {0, 0, 0, 0, "Simulation Options:", 3},
    {"number",    'n', "INT", 0, "Number of particles. Not in work when positions read from the file [default: 1]"},
    {"pbc",       301, 0, 0, "Use periodic boundary conditions [default: false]"},
    {"size",      401, "NUM", 0, "Box size [default: 10.0]"},
    {"depth",     501, "INT", 0, "SimpleTriangulation depth [default: 3]"},
    {"dt",        601, "NUM", 0, "Time step [default: 0.05]"},
    {"ttime",     602, "NUM", 0, "Total simulation time [default: 1.0]"},
    {"log-step",  603, "INT", 0, "[Log step interval [default: 1]"},
    {"int",       701, "STR", 0, "Integrator of equations of motion: velocity-verlet[vv], simple[se] [default: vv]"},

    {0,             0, 0, 0, "System Options:", 5},
    {0,           'a', "NUM", 0, "Repulsion parameter between bodies [default: 1.0]"},
    {"mass",      'm', "NUM", 0, "Mass of a particle [default: 1.0]"},
    {"gamma",     'k', "NUM", 0, "Spring constant [default: 1.0]"},
    {"mu",        801, "NUM", 0, "Viscosity coefficient [default: 1.0]"},
    {"dp",        802, "NUM", 0, "Osmotic pressure difference [default: 0.0]"},
    {"r-cut",     803, "NUM", 0, "Radius cut-off for pair interactions [default: 1.0]"},
    {0}
};


static int parse_opt (int key, char* arg, struct argp_state* state)
{
    /* Get the input argument from argp_parse, which
     * is a pointer to our arguments structure. */
    struct arguments* arguments = state->input;

    switch (key)
    {
        case ARGP_KEY_INIT:
            /* Default values. */
            arguments->silent = 0;
            arguments->verbose = 1;
            arguments->output_file = "./data/cells.xyz";
            arguments->input_file = "./data/cells.in";
            arguments->traj_file = "./data/render.py";
            arguments->log_file = "./data/biofilm.log";
            arguments->integrator_a = "vv";
            arguments->abort = 0;
            arguments->n_particles = 1;
            arguments->log_step = 1;
            arguments->r_cut = 1.0;
            arguments->dt = 0.05;
            arguments->ttime = 1.0;
            arguments->dp = 0.0;
            arguments->a = 1.0;
            arguments->d = 3;
            arguments->mass = 1.0;
            arguments->mu = 1.0;
            arguments->k = 1.0;
            arguments->L = 10.0;
            arguments->pbc = false;
            break;

        case 'q':
        case 's':
            arguments->silent = arg ? atoi (arg) : 1;
            arguments->verbose = 0;
            break;

        case 'v':
        case 'd':
            arguments->verbose = arg ? atoi (arg) : 1;
            arguments->silent = 0;
            break;

        case 'i':
            arguments->input_file = arg;
            break;

        case 'o':
            arguments->output_file = arg;
            break;

        case 't':
            arguments->traj_file = arg;
            break;

        case 'l':
            arguments->log_file = arg;
            break;

        case 'n':
            arguments->n_particles = arg ? atoi (arg) : 1;
            break;

        case 'a':
            arguments->a = arg ? strtod (arg, NULL) : 1.0;
            break;

        case 'm':
            arguments->mass = arg ? strtod (arg, NULL) : 1.0;
            break;
            
        case 'k':
            arguments->k = arg ? strtod (arg, NULL) : 1.0;
            break;

        case 301:
            arguments->pbc = true;
            break;
            
        case 401:
            arguments->L = arg ? strtod (arg, NULL) : 10.0;
            break;
            
        case 501:
            arguments->d = arg ? atoi (arg) : 3;
            break;

        case 601:
            arguments->dt = arg ? strtod (arg, NULL) : 0.001;
            break;
            
        case 602:
            arguments->ttime = arg ? strtod (arg, NULL) : 1.0;
            break;

        case 603:
            arguments->log_step = arg ? atoi (arg) : 1;
            break;
            
        case 701:
            arguments->integrator_a = arg;
            break;
            
         case 801:
            arguments->mu = arg ? strtod (arg, NULL) : 1.0;
            break;  
            
        case 802:
            arguments->dp = arg ? strtod (arg, NULL) : 0.0;
            break;
            
        case 803:
            arguments->r_cut = arg ?  strtod (arg, NULL) : 1.0;
            break;    

        case OPT_ABORT:
            arguments->abort = 1;
            break;

//      case ARGP_KEY_NO_ARGS:
//          argp_usage (state);

        case ARGP_KEY_ARG:
            arguments->strings = &state->argv[state->next];
            state->next = state->argc;
            break;

        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

Timer clocks[10];

int main(int argc, char** argv)
{   
    
    /* Initialize MT19937 Pseudo-random-number generator. */
    unsigned long init[4] = {0x123, 0x234, 0x345, 0x456}, length = 4;
    init_by_array(init, length);

    /* Parse our arguments; every option seen by parse_opt will
       be reflected in arguments. */
    struct arguments arguments;

    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    if (arguments.abort)
    {
        printf("PROGRAM FORCED TO >ABORT< \n");
        exit(123);
    }

    if (arguments.verbose && !arguments.silent)
    {
        printf ("OUTPUT_FILE = %s\n"
                "VERBOSE = %s\n"
                "SILENT = %s\n",
                arguments.output_file,
                arguments.verbose ? "yes" : "no",
                arguments.silent ? "yes" : "no");
    }
    
    ofstream os;
    os.open(arguments.log_file, ios::out | ios::trunc ); /* reset file */
    os.close();

    ofstream ost(arguments.traj_file, ios::out | ios::trunc);
    ost.close();
    
    int n_iter = (int) arguments.ttime / arguments.dt;
    int depth = arguments.d;
    

    clocks[0].tic();
    print_time();
    Vector3D va (0,0,0);
    Vector3D vb (1, 1,0);
    Vector3D vc (1,-1,0);
    Triangle t(va,vb,vc);
    /*-----------------------------------------------------------------------*/
    cout << "arguments.d " << arguments.d << endl;
    SimpleTriangulation sm(arguments.d);

    
    list<Triangle> tris = sm.triangulate();
    Cell cell(tris);
    double surf = cell.calcSurfaceArea();
    cout << "SURFACE AREA= " << surf <<  endl;
    cell.calcCM();
    double volume = cell.calcVolume();
    cout << "VOLUME= " << volume <<endl;
    
    cell.saveTriangulatedSurface(arguments.output_file);
    cell.saveRenderingScript(arguments.traj_file, arguments.output_file);
    
    cell.calcForces();
    
    clocks[0].toc();
    print_time();
    
    
    Vector3D v1(0, 0, 0);
    Vector3D v2(1, 1, 2);
    
    double gamma1 = 1.0;
    double gamma2 = 2.0;
    double R01 = 2.0;
    double R02 = sqrt(6.0);
    double R03 = 3.0;
    
    cout << HookeanForce::calcForce(v1,v2, R01,gamma1) <<endl;
    cout << HookeanForce::calcForce(v2,v1, R01,gamma1) <<endl;
    cout << HookeanForce::calcForce(v1,v2, R01,gamma2) <<endl;
    
    cout << HookeanForce::calcForce(v1,v2, R02,gamma1) <<endl;
    cout << HookeanForce::calcForce(v2,v1, R02,gamma1) <<endl;
    cout << HookeanForce::calcForce(v1,v2, R02,gamma2) <<endl;
    
    cout << HookeanForce::calcForce(v1,v2, R03,gamma1) <<endl;
    cout << HookeanForce::calcForce(v2,v1, R03,gamma1) <<endl;
    cout << HookeanForce::calcForce(v1,v2, R03,gamma2) <<endl;
    return 0;
}