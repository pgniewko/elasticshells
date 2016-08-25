#ifndef ARGUMENTS_H
#define	ARGUMENTS_H

struct arguments
{
    std::string render_file;            /* file arg to ‘--output’ */
    std::string traj_file;
    std::string box_file;
    std::string output_file;
    std::string surface_file;
    std::string stress_file;
    std::string topology_file;
    std::string lf_file;

    char** strings;               /* [string...] */

    int silent, verbose, debug, abort;   /* '-s', '-v', '-d','--abort' */

    char* files_prefix;
    char* output_dir;
    char* input_dir;
    std::string ob_config_file;
    std::string sch_config_file;

    char* integrator_a;
    char* tritype;
    char* model_type;

    int d;
    int log_step;
    int save_step;
    int box_step;
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
    double ttime;
    double r_vertex;                   /* vertex radius */
    double init_radius1;
    double init_radius2;
    double volume_scale;

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
    bool dynamics;
    bool nobending;
    bool const_volume;
    bool restart;
    int nb_flag;

    ulong seed;
};

#endif	/* ARGUMENTS_H */
