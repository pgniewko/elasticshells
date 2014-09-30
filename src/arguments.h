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
    char* surface_file;

    char* integrator_a;

    int d;
    int log_step;
    int save_step;
    int box_step;
    int vlist_step;
    int n_cells;
    int nsteps;

    double L;
    double a;
    double dt;
    double dp;
    double visc;
    double k;
    double mass;
    double ttime;
    double r_cut;                   /* cell1-cell2 */
    double r_bc;                    /* cell-box */
    double verlet_r;

    double bsx;
    double bsy;
    double bsz;
    double bsdx;
    double bsdy;
    double bsdz;
    //double bsxe;
    //double bsye;
    //double bsze;
    bool pbc;
    bool draw_box;
};

#endif	/* ARGUMENTS_H */