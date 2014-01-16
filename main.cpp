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
#include <error.h>     /* error */
#include <math.h>      /* log, sqrt */


#include "src/arguments.h"
#include "src/random.h"
#include "src/Simulator.h"
#include "src/Timer.h"

using namespace std;

const char* argp_program_version = "biofilm 0.0.9";
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
    {"wrap",     'w', "FILE",  0, "Coordinates wrapping mode:0 - image, 1 - real [default: 1]" },
    {"abort", OPT_ABORT, 0, 0, "Abort before showing any output"},
    
    {0, 0, 0, 0, "Simulation Options:", 3},
    {"int",       991,  "STR", 0, "Integrator of equations of motion: velocity-verlet[vv], trotter[trot] [default: trot]"},
    {"size",      444, "NUM", 0, "Box size [default: 10.0]"},
    {"n-iter",    666, "NUM", 0, "Number of time steps [default: 100]"},
    {"dt",        777, "NUM", 0, "Time step [default: 0.05]"},
    {"log-step",  888, "INT", 0, "[Log step interval [default: 1]"},
    {"number",    'n', "INT", 0, "Number of particles. Not in work when positions read from the file [default: 1]"},
    {"pbc",       301, 0, 0, "Use periodic boundary conditions [default: false]"},
    
    {0,             0, 0, 0, "System Options:", 5},
    { "r-cut",    555, "NUM", 0, "Radius cut-off for pair interactions [default: 1.0]"},
    {0,           'a', "NUM", 0, "Repulsion parameter between bodies [default: 25.0]"},
    {"gamma",     'g', "NUM", 0, "Medium viscosity [default: 4.5]"},
    {"sigma",     999, "NUM", 0, "Radius of a particle [default: 3.0]"},
    {"mass",      'm', "NUM", 0, "Mass of a particle [default: 1.0]"},
    {"pair-dist", 300, "NUM", 0, "[default: 1.0]"},
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
            arguments->output_file = "pos.coo";
            //arguments->input_file = "input.coo";
            arguments->input_file = "";
            arguments->traj_file = "traj.xyz";
            arguments->log_file = "log.txt";
            arguments->integrator_a = "trot";
            arguments->abort = 0;
            arguments->n_particles = 1;
            arguments->log_step = 1;
            arguments->L = 10.0;
            arguments->r_cut = 1.0;
            arguments->n_iter = 100;
            arguments->pair_dist = 1.0;
            arguments->dt = 0.05;
            arguments->a = 25.0;
            arguments->gamma = 4.5;
            arguments->sigma = 3.0;
            arguments->mass = 1.0;
            arguments->w = 1;
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
            
        case 'w':
            arguments->w = arg ? atoi (arg) : 1;
            break;
            
        case 'a':
            arguments->a = arg ? atoi (arg) : 25.0;
            break;
            
        case 'g':
            arguments->gamma = arg ? atoi (arg) : 4.5;
            break;

        case 'm':
            arguments->mass = arg ? atoi (arg) : 1.0;
            break;
            
        case 300:
            arguments->pair_dist = arg ? strtod (arg, NULL) : 1.0;
            break;
        
        case 301:
            arguments->pbc = true;
            break;
                 
        case 444:
            arguments->L = arg ? strtod (arg, NULL) : 10.0;
            break;
            
        case 777:
            arguments->dt = arg ? strtod (arg, NULL) : 0.001;
            break;
            
        case 666:
            arguments->n_iter = arg ?  atoi (arg) : 100;
            break;
            
        case 555:
            arguments->r_cut = arg ?  strtod (arg, NULL) : 1.0;
            break;    

        case 888:
            arguments->log_step = arg ? atoi (arg) : 1;
            break;

        case 991:
            arguments->integrator_a = arg;
            break;
            
        case 999:
            arguments->sigma = arg ? strtod (arg, NULL) : 3.0;
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
        error (123, 0, "ABORTED");
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

    //************ simulation is here
    clocks[0].tic();
    print_time();

    Box domain(arguments.L, arguments.L, arguments.L, Vector3D(0, 0, 0));

    Simulator model(arguments, domain);
    
    ofstream os;
    os.open(arguments.log_file, ios::out | ios::trunc ); /* reset file */
    os.close();
    
    ofstream ost(arguments.traj_file, ios::out | ios::trunc);
    ost.close();
    
    double sim_time = 0.0;

    for (int n = 0; n < arguments.n_iter; n++)
    {
        if ((n % arguments.log_step == 0) && arguments.verbose)
        {
            double T;
            Vector3D tot_P;
            model.state(T, tot_P);
            cout << "n=" << n << " t=" << sim_time << " temp=" << T << " total_momentum = " << tot_P.length() << endl;
            os.open(arguments.log_file, ios::app);
            os << sim_time << ' ' << T << '\n';
            os.close();
            
            model.write_pos_traj(arguments.traj_file, arguments.w);
        }

        clocks[2].tic();
        model.integrate();
        
        clocks[2].toc();
        sim_time += model.params.dt;
    }

    clocks[0].toc();
    int tt;

    for (tt = 0; tt < 3; tt++)
    {
        cout << " " << clocks[tt].time();
    }

    cout << '\n';
    
    model.write_pos(arguments.output_file, arguments.w);
    print_time();

    return 0;
}