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
    char* stress_file;
    char* ob_config_file;

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

    double E_cell;                         /* E cell-cell */
    double E_wall;                         /* E cell-wall */
    double nu;
    double thickness;
    double dt;
    double dp;
    double ddp;
    double eps;
    double visc;
//    double k;
    double mass;
    double ttime;
    double r_vertex;                   /* vertex radius */
    double verlet_r;
    double init_radius;
    double init_radius2;
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
    bool pbc;
    bool draw_box;
    bool osmotic_flag;
    bool scale_flag;
    int nb_flag;

    long seed;
};

#endif	/* ARGUMENTS_H */