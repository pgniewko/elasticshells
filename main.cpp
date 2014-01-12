/*
 * Author : Pawel Gniewek (UC Berkeley)
 * Email  : pawel.gniewek@berkeley.edu
 * License: BSD
 */

#include <stdio.h>     /* printf, fgets */
#include <argp.h>
#include <stdlib.h>     /* atoi,  strtod */
#include <string>
#include <iostream>
#include <error.h>     /* error */
#include <math.h>      /* log, sqrt */
#include <fstream>


#include "src/random.h"
#include "src/Simulator.h"
#include "src/Timer.h"

using namespace std;

struct Controls
{
    int transient;
    int update_step;
    char* fname;
};

Timer TT[10];


const char* argp_program_version = "biofilm 0.0.1";
const char* argp_program_bug_address = "<pawel.gniewek@berkeley.edu>";


static char doc[] =
    "General information about the software\
    \vThis part of the documentation comes *after* the options;";

static char args_doc[] = "[STRING...]";

struct arguments
{
    char** strings;               /* [string...] */
    int silent, verbose, abort;   /* ‘-s’, ‘-v’, ‘--abort’ */
    char* output_file;            /* file arg to ‘--output’ */
    int n;                        /* count arg to ‘--repeat’ */
    int t;
    double dt, mu;
    string cell_type;
};

#define OPT_ABORT  1            /* --abort */

static struct argp_option options[] =
{
    {0, 0, 0, 0, "Input/Output options:", 1},
    {"verbose",  'v', 0,       0, "Produce verbose output" },
    {"debug", 'd', 0, OPTION_ALIAS},
    {"quiet",    'q', 0,       0, "Don't produce any output" },
    {"silent",   's', 0,       OPTION_ALIAS },
    {"output",   'o', "FILE",  0, "Output to FILE instead of standard output" },
    {"abort",    OPT_ABORT, 0, 0, "Abort before showing any output"},
    {0, 0, 0, 0, "Simulation setup options:", -1},
    {"cell_type", 'c', "STR", 0, "[tumor, yeast, bacteria]"},
    { "time", 't', "NUM", OPTION_ARG_OPTIONAL, "..."},
    { "dt", 777, "NUM", OPTION_ARG_OPTIONAL, "..."},
    { "num", 'n', "INT", OPTION_ARG_OPTIONAL, "..."},
    //{ "num", 'n', "INT", OPTION_ARG_OPTIONAL, "..."}, <- when you must specify the value
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
            arguments->silent = 0;
            arguments->verbose = 0;
            arguments->output_file = "-";
            arguments->abort = 0;
            arguments->n = 1;
            break;

        case 'q':
        case 's':
            arguments->silent = 1;
            break;

        case 'v':
        case  'd':
            arguments->verbose = 1;
            break;

        case 'o':
            arguments->output_file = arg;
            break;

        case 'c':
            arguments->cell_type = arg;

        case 'n':
            arguments->n = arg ? atoi (arg) : 1;
            break;

        case 't':
            arguments->t = arg ? strtod (arg, NULL) : 1.0;
            break;

        case 777:
            arguments->dt = arg ? strtod (arg, NULL) : 0.001;

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
        error (10, 0, "ABORTED");
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
    TT[0].tic();
    print_time();
    Parameters Par;
    Controls C;
    // ../src/mydpd 0.05 100 input.coo out > log
    Par.Nparticles = 4000;
    Par.Nparticles = arguments.n;
    cout << Par.Nparticles << " " << arguments.n << endl;
    Par.L = 10;

    //Par.Niter=atoi(argv[2]);
    Par.Niter = 100;
    Par.rcut = 1.0;
    Par.pairdist = 1.0;
    //Par.tau=atof(argv[1]);
    Par.tau = 0.05;
    //Par.infname =  argv[3];
    Par.infname = "input.coo";
    Par.a = 25;
    Par.gamma = 4.5;
    Par.sigma = 3.0;
    Par.mass = 1.0;
    C.transient = 0;
    C.update_step = 1;
    //C.fname = argv[4];
    C.fname = "out";

    Distance domain(Par.L, Par.L, Par.L, Vector3D(0, 0, 0));

    Simulator model(Par, domain );

    ofstream os(C.fname);
    os.close();//reset file

    double t = 0.0;

    for (int n = 0; n < Par.Niter; n++)
    {
        if ((n % C.update_step == 0) && (n >= C.transient))
        {
            double T;
            Vector3D tot_P;
            model.state(T, tot_P);
            cout << "n=" << n << " t=" << t << " temp=" << T << " total_momentum = " << tot_P.length() << endl;
            ofstream os(C.fname, ios::app);
            os << t << ' ' << T << '\n';
            os.close();
        }

        TT[2].tic();
        model.integrate_trotter();
        TT[2].toc();
        t += model.Par.tau;
    }

    TT[0].toc();
    int tt;

    for (tt = 0; tt < 3; tt++)
    {
        cout << " " << TT[tt].time();
    }

    cout << '\n';
    model.write_pos("pos.coo");
    print_time();
    //*****************************************************************************

    return 0;
}
