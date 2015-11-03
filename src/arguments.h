#ifndef ARGUMENTS_H
#define	ARGUMENTS_H

//#include <string>
//using namespace std;

struct arguments
{
    std::string render_file;            /* file arg to ‘--output’ */
    std::string traj_file;
    std::string box_file;
    std::string output_file;
    std::string surface_file;
    std::string stress_file;

    char** strings;               /* [string...] */

    int silent, verbose, debug, abort;   /* '-s', '-v', '-d','--abort' */

    char* files_prefix;
    char* output_dir;
    char* input_dir;
    char* ob_config_file;
    char* sch_config_file;

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

    ulong seed;
};

#endif	/* ARGUMENTS_H */