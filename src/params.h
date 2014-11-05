#ifndef PARAMS_H
#define	PARAMS_H
//TODO - dokonczyc
struct params_t
{

    int log_step;
    int save_step;
    int box_step;
    int vlist_step;
    int n_cells;
    int nsteps;

    int d;
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


    bool draw_box;
    int nbhandler;
};

#endif	/* PARAMS_H */

