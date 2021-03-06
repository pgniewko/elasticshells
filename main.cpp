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
#include <climits>

#include "Environment.h"
#include "Timer.h"
#include "arguments.h"
#include "simulation/Simulator.h"
#include "utils/Logger.h"
#include "utils/LogManager.h"

utils::Logger shells_logs("elasticshells");

const char* argp_program_version = "elasticshells 0.9.9";
const char* argp_program_bug_address = "<gniewko.pablo@gmail.com>";

static char doc[] =
    "Copyright (C) 2014-2018, Pawel Gniewek. \nAll rights reserved.";

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

    {0,             0,       0, 0, "Simulation Options:", 4},
    {"number",    'n',   "INT", 0, "Init number of particles. Not in work when positions read from the file [default: 1]"},
    {"ns",        401,   "INT", 0, "Number of simulation steps [default: 10]"},
    {"dt",        402, "FLOAT", 0, "Integration time step [default: 0.001]"},
    {"pbc",       410,       0, 0, "Activate periodic boundary conditions [default: false]"},
    {"no-box",    411,       0, 0, "Deactivate box in rendering script - [default: true]"},
    {"tt",        412,   "STR", 0, "Triangulation type: Simple[simple], Platonic[plato], Random[rnd] [default: simple]"},
    {"depth",     413,   "INT", 0, "Triangulation depth [default: 5]"},
    {"plato",     414,   "INT", 0, "PlatonicTriangulation type: tetra[0], cube[1], octa[2], ico[3] [default: 0]"},
    {"bend",      417,       0, 0, "[default: false]"},
    {"restart",   419,       0, 0, "[default: false]"},
    {"analyze",   420,       0, 0, "[default: false]"},
    {"jam",       'j',       0, 0, "[default: false]"},

    {0,             0,       0, 0, "Material Properties:", 5},
    {"E-shell",   500, "FLOAT", 0, "Shell Young's modulus [default: 100.0]"},
    {"E-box",     501, "FLOAT", 0, "Box Young's modulus [default: 200.0]"},
    {"ir",        502, "FLOAT", 0, "Shells size at the initialization - lower limit [default:2.5"},
    {"ir2",       503, "FLOAT", 0, "Shells size at the initialization - upper limit [default:2.5]"},
    {"dp",        504, "FLOAT", 0, "Osmotic pressure [default: 1.0]"},
    {"osm",       505,       0, 0, "Volume dependent osmotic pressure [default:  false]"},
    {"rv",        506, "FLOAT", 0, "Radius of a single vertex [default: 0.25]"},
    {"ddp",       512, "FLOAT", 0, "Variation in osmotic pressure [UNIT=0.1 MPa] [default: 0.0]"},
    {"eps",       513, "FLOAT", 0, "Non osmotic volume fraction [default: 0.0]"},
    {"nu",        514, "FLOAT", 0, "Shell and box Poisson's ratio (the same for box and shell) [default: 0.5]"},
    {"th",        515, "FLOAT", 0, "Shell-wall thickness [UNIT=1 micron]  [default: 0.1]"},
    {"vol-f",     517,       0, 0, "Constant volume flag [default: false]"},
    {"min-force", 518, "FLOAT", 0, "Minimum force"},
    {"max-iter",  519,  "INT",  0, "Maximum number of iterations per minimization step"},

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

    struct arguments* arguments = static_cast<struct arguments*>(state->input);

    switch (key)
    {
        case ARGP_KEY_INIT:
            /* Default values. */
            arguments->silent = 0;
            arguments->verbose = 1;
            arguments->debug = 0;
            arguments->files_prefix = (char*)&"biofilm";
            arguments->output_dir = (char*)&"./output/";
            arguments->input_dir  = (char*)&"./input/";
            arguments->tritype = (char*)&"simple";
            arguments->d = 3;
            arguments->platotype = 0;
            arguments->n_shells = 1;
            arguments->nsteps = 10;
            arguments->E_shell = 100.0;
            arguments->E_wall = 200.0;
            arguments->nu = 0.5;
            arguments->thickness = 0.1;
            arguments->dt = 0.001;
            arguments->dp = 1.0;
            arguments->ddp = 0.0;
            arguments->eps = 0.0;
            arguments->r_vertex = 0.25;
            arguments->init_radius1 = 2.5;
            arguments->init_radius2 = 2.5;
            arguments->bsx = 10.0;
            arguments->bsy = 10.0;
            arguments->bsz = 10.0;
            arguments->bsdx = 0.0;
            arguments->bsdy = 0.0;
            arguments->bsdz = 0.0;
            arguments->bsxe = 5.0;
            arguments->bsye = 5.0;
            arguments->bsze = 5.0;
            arguments->pbc = false;
            arguments->osmotic_flag = false;
            arguments->bending = false;
            arguments->const_volume = false;
            arguments->restart = false;
            arguments->analyze = false;
            arguments->jam = false;
            arguments->seed = 0x123;
            arguments->min_force = 1e-6;
            arguments->max_iter = 10000;
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

        case 'n':
            arguments->n_shells = arg ? atoi (arg) : 1;
            break;

        case 401:
            arguments->nsteps = arg ? atoi (arg) : 100;
            break;

        case 402:
            arguments->dt = arg ? strtod (arg, NULL) : 0.001;
            break;

        case 410:
            arguments->pbc = true;
            break;

        case 412:
            arguments->tritype = arg;
            break;

        case 413:
            arguments->d = arg ? atoi (arg) : 5;
            break;

        case 414:
            arguments->platotype = arg ? atoi (arg) : 0;
            break;

        case 417:
            arguments->bending = true;
            break;

        case 419:
            arguments->restart = true;
            break;

        case 420:
            arguments->analyze = true;
            break;

        case 'j':
            arguments->jam = true;
            break;

        case 500:
            arguments->E_shell = arg ? strtod (arg, NULL) : 1500.0;
            break;

        case 501:
            arguments->E_wall = arg ? strtod (arg, NULL) : 2000.0;
            break;

        case 502:
            arguments->init_radius1 = arg ?  strtod (arg, NULL) : 2.5;
            break;

        case 503:
            arguments->init_radius2 = arg ?  strtod (arg, NULL) : 2.5;
            arguments->init_radius2 = std::max(arguments->init_radius1, arguments->init_radius2);
            break;

        case 504:
            arguments->dp = arg ? strtod (arg, NULL) : 0.0;
            break;

        case 505:
            arguments->osmotic_flag = true;
            break;

        case 506:
            arguments->r_vertex = arg ? strtod (arg, NULL) : 0.25;
            break;

        case 512:
            arguments->ddp = arg ? strtod (arg, NULL) : 0.0;
            break;

        case 513:
            arguments->eps = arg ? strtod (arg, NULL) : 0.0;
            break;

        case 514:
            arguments->nu = arg ? strtod (arg, NULL) : 0.5;
            break;

        case 515:
            arguments->thickness = arg ?  strtod (arg, NULL) : 0.1;
            break;

        case 517:
            arguments->const_volume = true;
            break;
            
        case 518:
            arguments->min_force = arg ? strtod (arg, NULL) : 1e-6;
            break;
        case 519:
            arguments->max_iter = arg ? atoi (arg) : 10000;
            break;

        case 601:
            arguments->bsx = arg ? strtod (arg, NULL) : 10.0;
            break;

        case 602:
            arguments->bsy = arg ? strtod (arg, NULL) : 10.0;
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
            
        case ARGP_KEY_ARG:
            arguments->strings = &state->argv[state->next];
            state->next = state->argc;
            break;

        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

void print_limits()
{
#ifdef DEBUG
    std::cout << "# =============================================" << std::endl;
    std::cout << "# Number of bits in byte:          " << CHAR_BIT << std::endl;
    std::cout << "# Maximum value of object(size_t): " << SIZE_MAX << std::endl;
    std::cout << "# Minimum value of int:            " << INT_MIN << std::endl;
    std::cout << "# Maximum value of int:            " << INT_MAX << std::endl;
    std::cout << "# Maximum value of signed int:     " << UINT_MAX << std::endl;
    std::cout << "# =============================================" << std::endl;
#endif
}
static struct argp argp = { options, parse_opt, args_doc, doc };

Timer clocks[10];
double simulation_time;

int main(int argc, char** argv)
{
    print_limits();

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
    
    arguments.traj_file    = std::string(arguments.output_dir) + std::string(arguments.files_prefix) + std::string(".xyz");
    arguments.box_file     = std::string(arguments.output_dir) + std::string(arguments.files_prefix) + std::string(".box.xyz");
    arguments.output_file  = std::string(arguments.output_dir) + std::string(arguments.files_prefix) + std::string(".out");
    arguments.topology_file = std::string(arguments.output_dir) + std::string(arguments.files_prefix) + std::string(".top");
    arguments.lf_file      = std::string(arguments.output_dir) + std::string(arguments.files_prefix) + std::string(".lf.xyz");

    arguments.ob_config_file = std::string(arguments.input_dir) + std::string("observers.config");
    arguments.sch_config_file = std::string(arguments.input_dir) + std::string("schedule.config");

    /* Initialize MT19937 Pseudo-random-number generator. */
    unsigned long init[4] = {arguments.seed, 0x234, 0x345, 0x456}, length = 4;
    init_by_array(init, length);
    shells_logs << utils::LogLevel::FILE << "TRAJ_FILE = "        << arguments.traj_file << "\n";
    shells_logs << utils::LogLevel::FILE << "BOX_FILE = "         << arguments.box_file << "\n";
    shells_logs << utils::LogLevel::FILE << "OUTPUT_FILE = "      << arguments.output_file << "\n";
    shells_logs << utils::LogLevel::FILE << "OBSERVERS_CONFIG = " << arguments.ob_config_file << "\n";

    if (arguments.n_shells == 0)
    {
        shells_logs << utils::LogLevel::INFO << "NUMBER OF SHELS IS ZERO (0). NOTHING TO DO !" << "\n";
    }
    else
    {
        clocks[0].tic();
        simulation_time = read_timer();
        Simulator simulator(arguments);

        if (arguments.analyze)
        {
            simulator.analyze();
        }
        else
        {
            if (arguments.restart)
            {
                simulator.restart();
            }
            else
            {
                simulator.init_shells(arguments.n_shells, arguments.init_radius1, arguments.init_radius2, arguments.jam);
            }

            simulator.simulate(arguments.nsteps);
        }

        clocks[0].toc();
        simulation_time = read_timer( ) - simulation_time;
    }

    shells_logs << utils::LogLevel::INFO << "TOTAL EXECUTION WALL-TIME = " << clocks[0].time() << " [s] \n";

    return (EXIT_SUCCESS);
}