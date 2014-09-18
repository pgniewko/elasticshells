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

#include "src/Timer.h"
//#include "src/random.h"
#include "src/arguments.h"
#include "src/geometry/Vector3D.h"
#include "src/simulation/Simulator.h"

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
    {"surf",     011, "FILE",  0, "Output to SURFACE-FILE instead of standard output [default: ... ]" },
    {"log",      'l', "FILE",  0, "Print log to FILE instead of standard output [default: ... ]" },
    {"xyz",      't', "FILE",  0, "Print trajectory to FILE [default: ... ]" },
    {"abort", OPT_ABORT, 0, 0, "Abort before showing any output"},

    {0, 0, 0, 0, "Simulation Options:", 3},
    {"number",    'n', "INT", 0, "Init number of particles. Not in work when positions read from the file [default: 1]"},
    {"pbc",       301, 0, 0, "Use periodic boundary conditions [default: false]"},
    {"dbox",      302, 0, 0, "Draw [default: true]"},
    {"size",      401, "NUM", 0, "Box size [default: 10.0]"},
    {"depth",     501, "INT", 0, "SimpleTriangulation depth [default: 3]"},
    {"dt",        601, "NUM", 0, "Time step [default: 0.01]"},
    {"ttime",     602, "NUM", 0, "Total simulation time [default: 1.0]"},
    {"log-step",  603, "INT", 0, "Log step interval [default: 10]"},
    {"ns",        604, "INT", 0, "Number of simulation steps [default: 100]"},
    {"save-step", 605, "INT", 0, "Save step interval [default: 10]"},
    {"box-step",  606, "INT", 0, "Box manipulation step interval [default: 10]"},
    {"vlist-step",607, "INT", 0, "Verlet-list step interval [default: 100]"},
    {"int",       701, "STR", 0, "Integrator of equations of motion: "
                                 "Forward-Euler[fe], Heun[hm], Runge-Kutta 2nd order[rk], Velocity-Verlet[vv] [default: fe]"},

    {0,             0, 0, 0, "System Options:", 5},
    {0,           'a', "NUM", 0, "Repulsion parameter between bodies [default: 1.0]"},
    {"mass",      'm', "NUM", 0, "Mass of a particle [default: 100.0]"},
    {"gamma",     'k', "NUM", 0, "Spring constant [default: 1.0]"},
    {"mu",        801, "NUM", 0, "Viscosity coefficient [default: 100.0]"},
    {"dp",        802, "NUM", 0, "Osmotic pressure [default: 0.0]"},
    {"r-cut",     803, "NUM", 0, "Radius cut-off for pair interactions [default: 1.0]"},
    {"bsx",       804, "NUM", 0, "X Box size [default: 10.0]"},
    {"bsy",       805, "NUM", 0, "Y Box size [default: 10.0]"},
    {"bsz",       806, "NUM", 0, "Z Box size [default: 10.0]"},
    {"verlet-r",  807, "NUM", 0, "Verlet radius [default: 2 times R_c]"},
    {"bsdx",      808, "NUM", 0, "dx of Box size [default: 0.0]"},
    {"bsdy",      809, "NUM", 0, "dy of Box size [default: 0.0]"},
    {"bsdz",      810, "NUM", 0, "dz of Box size [default: 0.0]"},
    {"rbc",       811, "NUM", 0, "Radius cut-off for cell-box interactions [default: 0.5]"},
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
            arguments->output_file = "./data/render.py";
            arguments->surface_file = "./data/surf.py";
            arguments->input_file = "./data/cells.in";
            arguments->traj_file = "./data/traj.xyz";
            arguments->log_file = "./data/biofilm.log";
            arguments->integrator_a = "vv";
            arguments->abort = 0;
            arguments->n_cells = 1;
            arguments->log_step = 10;
            arguments->save_step = 10;
            arguments->box_step = 10;
            arguments->vlist_step = 100;
            arguments->nsteps = 100;
            arguments->r_cut = 1.0;
            arguments->r_bc = 0.5;
            arguments->verlet_r = 2.0;
            arguments->dt = 0.01;
            arguments->ttime = 1.0;
            arguments->dp = 0.0;
            arguments->bsx = 10.0;
            arguments->bsy = 10.0;
            arguments->bsz = 10.0;
            arguments->bsdx = 0.0;
            arguments->bsdy = 0.0;
            arguments->bsdz = 0.0;
            arguments->a = 1.0;
            arguments->d = 3;
            arguments->mass = 100.0;
            arguments->visc = 100.0;
            arguments->k = 1.0;
            arguments->L = 10.0;
            arguments->pbc = false;
            arguments->draw_box = true;
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
            arguments->n_cells = arg ? atoi (arg) : 1;
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
            
        case 011:
            arguments->surface_file = arg;
            break;    

        case 301:
            arguments->pbc = true;
            break;
            
        case 302:
            arguments->draw_box = true;
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
            
        case 604:
            arguments->nsteps = arg ? atoi (arg) : 100;
            break;
            
        case 605:
            arguments->save_step = arg ? atoi (arg) : 1;
            break;
            
        case 606:
            arguments->box_step = arg ? atoi (arg) : 1;
            break;
            
        case 607:
            arguments->vlist_step = arg ? atoi (arg) : 100;
            break;    
            
        case 701:
            arguments->integrator_a = arg;
            break;
            
         case 801:
            arguments->visc = arg ? strtod (arg, NULL) : 1.0;
            break;  
            
        case 802:
            arguments->dp = arg ? strtod (arg, NULL) : 0.0;
            break;
            
        case 803:
            arguments->r_cut = arg ?  strtod (arg, NULL) : 1.0;
            break;
           
        case 804:
            arguments->bsx = arg ?  strtod (arg, NULL) : 10.0;
            break;
            
        case 805:
            arguments->bsy = arg ?  strtod (arg, NULL) : 10.0;
            break;
            
        case 806:
            arguments->bsz = arg ?  strtod (arg, NULL) : 10.0;
            break;
            
        case 807:
            arguments->verlet_r = arg ?  strtod (arg, NULL) : 2.0;
            break;
            
        case 808:
            arguments->bsdx = arg ?  strtod (arg, NULL) : 0.0;
            break;
            
        case 809:
            arguments->bsdy = arg ?  strtod (arg, NULL) : 0.0;
            break;
            
        case 810:
            arguments->bsdz = arg ?  strtod (arg, NULL) : 0.0;
            break;       
            
        case 811:
            arguments->r_bc = arg ?  strtod (arg, NULL) : 0.5;
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
    
//    /* Initialize MT19937 Pseudo-random-number generator. */
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
    
    clocks[0].tic();
    print_time();
    
    

    Simulator simulator(arguments);
    simulator.initCells(arguments.n_cells, 1.5);
    //simulator.addCell(1.5);
    //simulator.addCell(1.5);
    
    
    //cout << max(122,123) << endl;
    
    //Vector3D vel(-.5,-.5,-.5);
    //Vector3D vel(-.0,-.0,-.0);
    //Vector3D shift(4,4,4);
    //Vector3D shift(0,-3.4,0);
    //simulator.addCellVel(-vel, 0);
    //simulator.addCellVel(vel, 1);
    //simulator.moveCell(shift, 1);
    //simulator.calcForces();
    simulator.simulate(arguments.nsteps);
    //simulator.printCell(0);
    
    
    //cout << "cell #1 mass = " << simulator.getCell(0).getMass() << endl;
    //cout << "cell #2 mass = " << simulator.getCell(1).getMass() << endl;
    
    return (EXIT_SUCCESS);
}