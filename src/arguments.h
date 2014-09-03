#ifndef ARGUMENTS_H
#define	ARGUMENTS_H

//#include <string>
using namespace std;

struct arguments
{
    char** strings;               /* [string...] */

    int silent, verbose, abort;   /* ‘-s’, ‘-v’, ‘--abort’ */

    char* input_file;             /* file arg to ‘--input’ */
    char* output_file;            /* file arg to ‘--output’ */
    char* traj_file;              /* file arg to ‘--traj’ */
    char* log_file;               /* file arg to ‘--log’ */

    char* integrator_a;

    int d;
    int log_step;
    int n_particles;
    int nsteps;

    double L;
    double a;
    double dt;
    double dp;
    double mu;
    double k;
    double mass;
    double ttime;
    double r_cut;
    double bs;

    bool pbc;
    bool draw_box;
};

#endif	/* ARGUMENTS_H */