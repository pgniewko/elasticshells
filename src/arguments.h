#ifndef ARGUMENTS_H
#define	ARGUMENTS_H

//#include <string>
using namespace std;

struct arguments
{
    char** strings;               /* [string...] */
    
    
    int silent, verbose, abort;   /* ‘-s’, ‘-v’, ‘--abort’ */
    
    
    char* input_file;          /* file arg to ‘--input’ */
    char* output_file;         /* file arg to ‘--output’ */
    char* traj_file;         /* file arg to ‘--traj’ */
    char* log_file;            /* file arg to ‘--log’ */
    
    char* integrator_a;
    
    int w;
    int log_step;
    int n_particles;
    int n_iter;
    
    double L;
    double a;
    double dt;
    double mu;
    double gamma;
    double sigma;
    double mass;
    double r_cut;
    double pair_dist;
    
    bool pbc;
};

#endif	/* ARGUMENTS_H */