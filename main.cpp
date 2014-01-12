/*
 * Author : Pawel Gniewek (UC Berkeley)
 * Email  : pawel.gniewek@berkeley.edu
 * License: BSD
 */

#include <stdio.h>     /* printf, fgets */
#include <argp.h>
#include <stdlib.h>     /* atoi,  strtod */
//#include <string>
#include <iostream>
#include <error.h>     /* error */
#include <math.h>      /* log, sqrt */
#include <fstream>


#include "src/arguments.h"
#include "src/random.h"
#include "src/Simulator.h"
#include "src/Timer.h"

using namespace std;

const char* argp_program_version = "biofilm 0.0.1";
const char* argp_program_bug_address = "<pawel.gniewek@berkeley.edu>";


static char doc[] =
    "General information about the software\
    \vThis part of the documentation comes *after* the options;";

static char args_doc[] = "[STRING...]";

#define OPT_ABORT  1            /* --abort */

static struct argp_option options[] =
{
    {0, 0, 0, 0, "Input/Output options:", 1},
    {"verbose",  'v', 0, 0, "Produce verbose output" },
    {"debug", 'd', 0, OPTION_ALIAS},
    {"quiet",    'q', 0, 0, "Don't produce any output" },
    {"silent",   's', 0, OPTION_ALIAS},
    {"input",   'i', "FILE",  0, "Input from FILE" },
    {"output",   'o', "FILE",  0, "Output to FILE instead of standard output" },
    {"log",   'l', "FILE",  0, "Print log to FILE instead of standard output" },
    {"abort", OPT_ABORT, 0, 0, "Abort before showing any output"},
    
    {0, 0, 0, 0, "Simulation setup options:", -1},
    {"cell-type", 'c', "STR", 0, "[tumor, yeast, bacteria]"},
    { "time", 't', "NUM", 0, "..."},
    { "dl", 444, "NUM", 0, "..."},
    { "n-iter", 666, "NUM", 0, "..."},
    { "dt", 777, "NUM", 0, "..."},
    { "step", 888, "INT", 0, "..."},
    { "number", 'n', "INT", 0, "..."},
    {"pbc", 111, 0, 0, "Periodic boundary conditions"},
    {0, 0, 0, 0, "Simulation run options:", -1},
    { "r-cut", 555, "NUM", 0, "..."},
    {0, 'a', "NUM", 0, "..."},
    {"gamma", 'g', "NUM", 0, "..."},
    {"sigma", 999, "NUM", 0, "..."},
    {"mass", 'm', "NUM", 0, "..."},
    {"pair-dist", 011, "NUM", 0, "..."},
    //{ "n-iter", 666, "NUM", OPTION_ARG_OPTIONAL, "..."},
    {0}
};


static int parse_opt (int key, char* arg, struct argp_state* state)
{
    /* Get the input argument from argp_parse, which we
       know is a pointer to our arguments structure. */
    struct arguments* arguments = state->input;

    switch (key)
    {
        case ARGP_KEY_INIT:
            /* Default values. */
            arguments->silent = 1;
            arguments->verbose = 0;
            arguments->output_file = "pos.coo";
            arguments->input_file = "input.coo";
            arguments->log_file = "log.txt";
            arguments->abort = 0;
            arguments->n_particles = 1;
            arguments->update_step = 1;
            arguments->L = 10.0;
            arguments->rcut = 1.0;
            arguments->n_iter = 100;
            arguments->pair_dist = 1.0;
            arguments->dt = 0.05;
            arguments->a = 25;
            arguments->gamma = 4.5;
            arguments->sigma = 3.0;
            arguments->mass = 1.0;
            arguments->pbc = false;
            break;

        case 'q':
        case 's':
            arguments->silent = 1;
            break;

        case 'v':
        case 'd':
            arguments->verbose = 1;
            break;

        case 'i':
            arguments->input_file = arg;
            break;
            
        case 'o':
            arguments->output_file = arg;
            break;
            
        case 'l':
            arguments->log_file = arg;
            break;
            
        case 'c':
            arguments->cell_type = arg;

        case 'n':
            arguments->n_particles = arg ? atoi (arg) : 1;
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
            
        case 011:
            arguments->pair_dist = arg ? strtod (arg, NULL) : 1.0;
            break;
        
        case 111:
            arguments->pbc = false;
            break;
                 
        case 777:
            arguments->dt = arg ? strtod (arg, NULL) : 0.001;
            break;
            
        case 666:
            arguments->n_iter = arg ?  atoi (arg) : 100;

        case 888:
            arguments->update_step = arg ? atoi (arg) : 1;
            
        case 999:
            arguments->sigma = arg ? strtod (arg, NULL) : 3.0;
            
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

    //************simulatoin is here
    clocks[0].tic();
    print_time();

    Distance domain(arguments.L, arguments.L, arguments.L, Vector3D(0, 0, 0));

    Simulator model(arguments, domain );
    
    ofstream os(arguments.output_file);
    os.close();//reset file
    double t = 0.0;

    for (int n = 0; n < arguments.n_iter; n++)
    {
        if ((n % arguments.update_step == 0) && arguments.verbose == 0)
        {
            double T;
            Vector3D tot_P;
            model.state(T, tot_P);
            cout << "n=" << n << " t=" << t << " temp=" << T << " total_momentum = " << tot_P.length() << endl;
            ofstream os(arguments.output_file, ios::app);
            os << t << ' ' << T << '\n';
            os.close();
        }

        clocks[2].tic();
        model.integrate_trotter();
        clocks[2].toc();
        t += model.params.dt;
    }

    clocks[0].toc();
    int tt;

    for (tt = 0; tt < 3; tt++)
    {
        cout << " " << clocks[tt].time();
    }

    cout << '\n';
    
    model.write_pos(arguments.output_file);
    print_time();

    return 0;
}
