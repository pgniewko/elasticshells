#ifndef ARGUMENTS_H
#define	ARGUMENTS_H

struct arguments
{
    std::string traj_file;
    std::string box_file;
    std::string output_file;
    std::string topology_file;
    std::string lf_file;

    char** strings;               /* [string...] */

    int silent, verbose, debug; //, abort;   /* '-s', '-v', '-d','--abort' */

    char* files_prefix;
    char* output_dir;
    char* input_dir;
    std::string ob_config_file;
    std::string sch_config_file;

    char* tritype;

    int d;
    int n_shells;
    int nsteps;
    int platotype;

    double E_shell;                         /* E cell-cell */
    double E_wall;                         /* E cell-wall */
    double nu;
    double thickness;
    double dt;
    double dp;
    double ddp;
    double eps;
    double r_vertex;                   /* vertex radius */
    double init_radius1;
    double init_radius2;

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
    bool osmotic_flag;
    bool bending;
    bool const_volume;
    bool restart;
    bool analyze;
    bool jam;

    unsigned long seed;
    double min_force;
    int max_iter;
    
    bool ellipsoid;
    double rx;
    double ry;
    double rz;
    int n_verts;
};

#endif	/* ARGUMENTS_H */
