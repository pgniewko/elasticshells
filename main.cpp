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
#include <string>

#include "Environment.h"
#include "src/Timer.h"
#include "src/arguments.h"
#include "src/geometry/Vector3D.h"
#include "src/simulation/Simulator.h"

#include "utils/Logger.h"
#include "src/utils/LogManager.h"

utils::Logger biofilm_logs("biofilm");

const char* argp_program_version = "biofilm 0.5";
const char* argp_program_bug_address = "<pawel.gniewek@berkeley.edu>";

static char doc[] =
    "General information about the software goes *here*\
    \vCopyright (C) 2014-2015, Pawel Gniewek. \nAll rights reserved.";

static char args_doc[] = "";

#define OPT_ABORT  1            /* --abort ^A*/

static struct argp_option options[] =
{
    {0,               0,      0,  0, "I/O Options:", 3},
    {"verbose",     'v',  "INT",  OPTION_ARG_OPTIONAL, "Produce verbose output for values [default: 1]" },
    {"debug",       'd',      0,  OPTION_ALIAS},
    {"quiet",       'q',  "INT",  OPTION_ARG_OPTIONAL, "Don't produce any output [default: 1]" },
    {"silent",      's',      0,  OPTION_ALIAS},
    {"input",       'i', "FILE",  0, "Input from FILE [default: ...]" },
    {"out-dir",     301, "FILE",  0, "... [default: ./output]" },
    {"in-dir",      302,  "STR",  0, "... [default: ./input]" },
    {"prefix",      303,  "STR",  0, "... [default: biofilm]" },
    {"seed",        304, "LONG",  0, "Random generator seed [default: 0x123] " },
    {"oc",          305, "FILE",  0, "Input file to configure observers [default: ./observers.config ]" },
    {"sch",         306, "FILE",  0, "Input file to configure compression schedule [default: ./schedule.config ]" },
    {"abort", OPT_ABORT, 0, 0, "Abort before showing any output"},

    {0,             0,       0, 0, "Simulation Options:", 4},
    {"number",    'n',   "INT", 0, "Init number of particles. Not in work when positions read from the file [default: 1]"},
    {"time",      't', "FLOAT", 0, "Total simulation time [default: 1.0]"},
    {"ns",        401,   "INT", 0, "Number of simulation steps [default: 10]"},
    {"dt",        402, "FLOAT", 0, "Time step [default: 0.001]"},
    {"int",       403,   "STR", 0, "Integrator of equations of motion: Forward-Euler[fe], Heun[hm], Runge-Kutta 2nd order[rk] [default: fe]"},
    {"nb",        404,   "INT", 0, "Nb interaction handler: Naive O(N^2)[0], Verlet-list[1], Linked-domains[2] [default: 0]"},
    {"log-step",  405,   "INT", 0, "Log step interval [default: 10]"},
    {"save-step", 406,   "INT", 0, "Save step interval [default: 1]"},
    {"box-step",  407,   "INT", 0, "Box manipulation step interval [default: 10]"},
    {"verlet-r",  409, "FLOAT", 0, "Verlet radius times r_vertex [default: 3]"},
    {"pbc",       410,       0, 0, "Activate periodic boundary conditions [default: false]"},
    {"no-box",    411,       0, 0, "Deactivate box in rendering script - [default: true]"},
    {"tt",        412,   "STR", 0, "Triangulation type: Simple[simple], Platonic[plato] [default: simple]"},
    {"depth",     413,   "INT", 0, "Triangulation depth [default: 3]"},
    {"plato",     414,   "INT", 0, "PlatonicTriangulation type: tetra[0], cube[1], ico[1], oct[2] [default: 0]"},
    {"scale",     415,       0, 0, "Scale the saved coordinates upon compression [default: false]"},

    {0,             0,       0, 0, "Cell Options:", 5},
    {"ecc",       500, "FLOAT", 0, "Cell-wall Young's modulus [UNIT=0.1 MPa] [default: 1500.0]"},
    {"ecw",       501, "FLOAT", 0, "Box Young's modulus [UNIT=0.1 MPa] [default: 2000.0]"},
    {"ir",        502, "FLOAT", 0, "Cells size at the initialization - lower limit [default:2.5"},
    {"mu",        503, "FLOAT", 0, "Viscosity coefficient [default: 100.0]"},
    {"dp",        504, "FLOAT", 0, "Osmotic pressure [default: 0.0]"},
    {"osm",       505,       0, 0, "Volume dependent osmotic pressure [default:  false]"},
    {"rv",        506, "FLOAT", 0, "Radius of a single vertex [default: 0.25]"},
    {"gr",        507, "FLOAT", 0, "Cell growth rate [default: 0.0]"},
    {"vc",        508, "FLOAT", 0, "Minimum volume for budding[default: 60.0]"},
    {"bd",        509, "FLOAT", 0, "Bud-neck diameter [default: 0.5]"},
    {"dr",        510, "FLOAT", 0, "Bud volume and Mother cell volume ratio at division [default: 0.7]"},
    {"ir2",       511, "FLOAT", 0, "Cells size at the initialization - upper limit [default:2.5]"},
    {"ddp",       512, "FLOAT", 0, "Variation in osmotic pressure [UNIT=0.1 MPa] [default: 0.0]"},
    {"eps",       513, "FLOAT", 0, "Non osmotic volume fraction [default: 0.0]"},
    {"nu",        514, "FLOAT", 0, "Cell and box Poisson's ratio (the same for box and cell) [default: 0.5]"},
    {"th",        515, "FLOAT", 0, "Cell-wall thickness [UNIT=1 micron]  [default: 0.1]"},

    {0,             0,       0, 0, "Box options:", 6},
    {"bsx",       601, "FLOAT", 0, "X Box size [default: 10.0]"},
    {"bsy",       602, "FLOAT", 0, "Y Box size [default: 10.0]"},
    {"bsz",       603, "FLOAT", 0, "Z Box size [default: 10.0]"},
    {"bsdx",      604, "FLOAT", 0, "dx of Box size [default: 0.0]"},
    {"bsdy",      605, "FLOAT", 0, "dy of Box size [default: 0.0]"},
    {"bsdz",      606, "FLOAT", 0, "dz of Box size [default: 0.0]"},
    {"bsxe",      607, "FLOAT", 0, "X end of Box size [default: 10.0]"},
    {"bsye",      608, "FLOAT", 0, "Y end of Box size [default: 10.0]"},
    {"bsze",      609, "FLOAT", 0, "Z end of Box size [default: 10.0]"},
    {0}
};


static int parse_opt (int key, char* arg, struct argp_state* state)
{
    /* Get the input argument from argp_parse, which
     * is a pointer to our arguments structure. */

    //struct arguments* arguments = (struct arguments*) state->input;
    struct arguments* arguments = static_cast<struct arguments*>(state->input);

    switch (key)
    {
        case ARGP_KEY_INIT:
            /* Default values. */
            arguments->silent = 0;
            arguments->verbose = 1;
            arguments->debug = 0;
            arguments->abort = 0;
            arguments->files_prefix = (char*)&"biofilm";
            arguments->output_dir = (char*)&"./output/";
            arguments->input_dir  = (char*)&"./input/";
            arguments->ob_config_file = (char*)&"./input/observers.config";
            arguments->sch_config_file = (char*)&"./input/schedule.config";
            arguments->integrator_a = (char*)&"fe";
            arguments->tritype = (char*)&"simple";
            arguments->d = 3;
            arguments->platotype = 0;
            arguments->log_step = 10;
            arguments->save_step = 1;
            arguments->box_step = 10;
            arguments->n_cells = 1;
            arguments->nsteps = 10;
            arguments->E_cell = 1500.0;
            arguments->E_wall = 2000.0;
            arguments->nu = 0.5;
            arguments->thickness = 0.1;
            arguments->dt = 0.001;
            arguments->dp = 0.0;
            arguments->ddp = 0.0;
            arguments->eps = 0.0;
            arguments->visc = 100.0;
            arguments->ttime = 1.0;
            arguments->r_vertex = 0.25;
            arguments->verlet_f = 3.0;
            arguments->init_radius = 2.5;
            arguments->init_radius2 = 2.5;
            arguments->growth_rate = 0.0;
            arguments->vc = 60.0;
            arguments->bud_d = 0.5;
            arguments->div_ratio = 0.7;
            arguments->bsx = 10.0;
            arguments->bsy = 10.0;
            arguments->bsz = 10.0;
            arguments->bsdx = 0.0;
            arguments->bsdy = 0.0;
            arguments->bsdz = 0.0;
            arguments->bsxe = 10.0;
            arguments->bsye = 10.0;
            arguments->bsze = 10.0;
            arguments->pbc = false;
            arguments->draw_box = true;
            arguments->osmotic_flag = false;
            arguments->scale_flag = false;
            arguments->nb_flag = 0;
            arguments->seed = 0x123;
            break;

        case 'q':
        case 's':
            arguments->silent = arg ? atoi (arg) : 1;
            arguments->verbose = 0;
            break;

        case 'v':
            arguments->verbose = arg ? atoi (arg) : 1;
            arguments->silent = 0;
            break;

        case 'd':
            arguments->debug = arg ? atoi (arg) : 1;
            arguments->silent = 0;
            break;

        case 'i':
            break;

        case 301:
            arguments->output_dir =  arg;
            break;

        case 302:
            arguments->input_dir  =  arg;
            break;

        case 303:
            arguments->files_prefix =  arg;
            break;

        case 304:
            arguments->seed =  arg ? atol (arg) : 0x123;
            break;

        case 305:
            arguments->ob_config_file = arg;
            break;

        case 306:
            arguments->sch_config_file = arg;
            break;

        case 'n':
            arguments->n_cells = arg ? atoi (arg) : 1;
            break;

        case 't':
            arguments->ttime = arg ? strtod (arg, NULL) : 1.0;
            break;

        case 401:
            arguments->nsteps = arg ? atoi (arg) : 100;
            break;

        case 402:
            arguments->dt = arg ? strtod (arg, NULL) : 0.001;
            break;

        case 403:
            arguments->integrator_a = arg;
            break;

        case 404:
            arguments->nb_flag = arg ? atoi (arg) : 0;
            break;

        case 405:
            arguments->log_step = arg ? atoi (arg) : 10;
            break;

        case 406:
            arguments->save_step = arg ? atoi (arg) : 1;
            break;

        case 407:
            arguments->box_step = arg ? atoi (arg) : 10;
            break;

        case 409:
            arguments->verlet_f = arg ?  strtod (arg, NULL) : 3.0;
            arguments->verlet_f = std::max(arguments->verlet_f, 1.0);
            break;

        case 410:
            arguments->pbc = true;
            break;

        case 411:
            arguments->draw_box = false;
            break;

        case 412:
            arguments->tritype = arg;
            break;

        case 413:
            arguments->d = arg ? atoi (arg) : 3;
            break;

        case 414:
            arguments->platotype = arg ? atoi (arg) : 0;
            break;

        case 415:
            arguments->scale_flag = true;
            break;

        case 500:
            arguments->E_cell = arg ? strtod (arg, NULL) : 1500.0;
            break;

        case 501:
            arguments->E_wall = arg ? strtod (arg, NULL) : 2000.0;
            break;

        case 502:
            arguments->init_radius = arg ?  strtod (arg, NULL) : 2.5;
            break;

        case 503:
            arguments->visc = arg ? strtod (arg, NULL) : 100.0;
            break;

        case 504:
            arguments->dp = arg ? strtod (arg, NULL) : 0.0;
            break;

        case 505:
            arguments->osmotic_flag = true;
            break;

        case 506:
            arguments->r_vertex = arg ?  strtod (arg, NULL) : 0.25;
            break;

        case 507:
            arguments->growth_rate = arg ?  strtod (arg, NULL) : 0.0;
            break;

        case 508:
            arguments->vc = arg ?  strtod (arg, NULL) : 60.0;
            break;

        case 509:
            arguments->bud_d = arg ?  strtod (arg, NULL) : 0.5;
            break;

        case 510:
            arguments->div_ratio = arg ?  strtod (arg, NULL) : 0.7;
            break;

        case 511:
            arguments->init_radius2 = arg ?  strtod (arg, NULL) : 2.5;
            arguments->init_radius2 = std::max(arguments->init_radius, arguments->init_radius2);
            break;

        case 512:
            arguments->ddp = arg ?  strtod (arg, NULL) : 0.0;
            break;

        case 513:
            arguments->eps = arg ?  strtod (arg, NULL) : 0.0;
            break;

        case 514:
            arguments->nu = arg ?  strtod (arg, NULL) : 0.5;
            break;

        case 515:
            arguments->thickness = arg ?  strtod (arg, NULL) : 0.1;
            break;

        case 601:
            arguments->bsx = arg ?  strtod (arg, NULL) : 10.0;
            break;

        case 602:
            arguments->bsy = arg ?  strtod (arg, NULL) : 10.0;
            break;

        case 603:
            arguments->bsz = arg ?  strtod (arg, NULL) : 10.0;
            break;

        case 604:
            arguments->bsdx = arg ?  strtod (arg, NULL) : 0.0;
            break;

        case 605:
            arguments->bsdy = arg ?  strtod (arg, NULL) : 0.0;
            break;

        case 606:
            arguments->bsdz = arg ?  strtod (arg, NULL) : 0.0;
            break;

        case 607:
            arguments->bsxe = arg ?  strtod (arg, NULL) : 10.0;
            break;

        case 608:
            arguments->bsye = arg ?  strtod (arg, NULL) : 10.0;
            break;

        case 609:
            arguments->bsze = arg ?  strtod (arg, NULL) : 10.0;
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
double simulation_time;

int main(int argc, char** argv)
{
    print_time();

    if ( argc <= 1 )
    {
        argp_help(&argp, stdout, ARGP_HELP_SEE, argv[0]);
        exit(EXIT_FAILURE);
    }

    /* Parse our arguments; every option seen by parse_opt will
       be reflected in arguments. */
    struct arguments arguments;
    argp_parse (&argp, argc, argv, 0, 0, &arguments);

    if (arguments.debug && !arguments.silent)
    {
        utils::LogManager::set_level("FINEST");
    }
    else if (arguments.verbose && !arguments.silent)
    {
        utils::LogManager::set_level("INFO");
    }
    else if (arguments.silent)
    {
        utils::LogManager::set_level("SEVERE");
    }

    if (arguments.abort)
    {
        biofilm_logs << utils::LogLevel::SEVERE << "PROGRAM FORCED TO *ABORT*\n";
        exit(EXIT_FAILURE);
    }

    arguments.render_file  = std::string(arguments.output_dir) + std::string(arguments.files_prefix) + std::string(".py");
    arguments.traj_file    = std::string(arguments.output_dir) + std::string(arguments.files_prefix) + std::string(".xyz");
    arguments.box_file     = std::string(arguments.output_dir) + std::string(arguments.files_prefix) + std::string(".box.xyz");
    arguments.output_file  = std::string(arguments.output_dir) + std::string(arguments.files_prefix) + std::string(".out");
    arguments.surface_file = std::string(arguments.output_dir) + std::string(arguments.files_prefix) + std::string(".surface.py");
    arguments.stress_file  = std::string(arguments.output_dir) + std::string(arguments.files_prefix) + std::string(".stress.py");

    /* Initialize MT19937 Pseudo-random-number generator. */
    unsigned long init[4] = {arguments.seed, 0x234, 0x345, 0x456}, length = 4;
    init_by_array(init, length);
    biofilm_logs << utils::LogLevel::FILE << "RENDER_FILE = "      << arguments.render_file << "\n";
    biofilm_logs << utils::LogLevel::FILE << "TRAJ_FILE = "        << arguments.traj_file << "\n";
    biofilm_logs << utils::LogLevel::FILE << "BOX_FILE = "         << arguments.box_file << "\n";
    biofilm_logs << utils::LogLevel::FILE << "OUTPUT_FILE = "      << arguments.output_file << "\n";
    biofilm_logs << utils::LogLevel::FILE << "SURFACE_FILE = "     << arguments.surface_file << "\n";
    biofilm_logs << utils::LogLevel::FILE << "STRESS_FILE = "      << arguments.stress_file << "\n";
    biofilm_logs << utils::LogLevel::FILE << "OBSERVERS_CONFIG = " << arguments.ob_config_file << "\n";

    clocks[0].tic();
    simulation_time = read_timer();
    Simulator simulator(arguments);
    simulator.initCells(arguments.n_cells, arguments.init_radius, arguments.init_radius2);
    simulator.simulate(arguments.nsteps);
    clocks[0].toc();
    simulation_time = read_timer( ) - simulation_time;

#ifdef _OPENMP
    int gt = omp_get_max_threads();
    int ncpu = omp_get_num_procs();
    biofilm_logs << utils::LogLevel::INFO << "TOTAL EXECUTION WALL-TIME = " << simulation_time << " [s] \n";
    biofilm_logs << utils::LogLevel::INFO << "TOTAL EXECUTION CPU-TIME[" << gt << "(n_th)/" << ncpu << "(n_pcu)] = " << clocks[0].time() << " [s] \n";
#else
    biofilm_logs << utils::LogLevel::INFO << "TOTAL EXECUTION WALL-TIME = " << clocks[0].time() << " [s] \n";
#endif


    return (EXIT_SUCCESS);
}