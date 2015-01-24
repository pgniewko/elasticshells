#ifndef ARGUMENTS_H
#define	ARGUMENTS_H

//using namespace std;

struct arguments
{
    char** strings;               /* [string...] */

    int silent, verbose, debug, abort;   /* '-s', '-v', '-d','--abort' */

    char* render_file;            /* file arg to ‘--output’ */
    char* traj_file;              /* file arg to ‘--traj’ */
    char* output_file;            /* file arg to ‘--log’ */
    char* surface_file;

    char* integrator_a;
    char* tritype;

    int d;
    int log_step;
    int save_step;
    int box_step;
    int vlist_step;
    int n_cells;
    int nsteps;
    int platotype;

    double ecc;                         /* E cell-cell */
    double dt;
    double dp;
    double visc;
    double k;
    double mass;
    double ttime;
    double r_vertex;                   /* vertex radius */
    double verlet_r;
    double init_radius;
    double growth_rate;
    double vc;
    double bud_d;
    double div_ratio;

    double bsx;
    double bsy;
    double bsz;
    double bsdx;
    double bsdy;
    double bsdz;
    double bsxe;
    double bsye;
    double bsze;
    double ecw;                         /* E cell-wall */
    bool pbc;
    bool draw_box;
    bool osmotic_flag;
    int nb_flag;
};

#endif	/* ARGUMENTS_H */