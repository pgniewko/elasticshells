#ifndef ARGUMENTS_H
#define	ARGUMENTS_H

#include <string>
using namespace std;

struct arguments
{
    char** strings;               /* [string...] */
    
    
    int silent, verbose, abort;   /* ‘-s’, ‘-v’, ‘--abort’ */
    
    
    char* input_file;            /* file arg to ‘--output’ */
    char* output_file;            /* file arg to ‘--output’ */
    char* log_file;            /* file arg to ‘--output’ */
    
    
    int update_step;
    double mu;
    string cell_type;
    int n_particles;
    double L;
    int n_iter;
    double rcut;
    double dt;
    double pair_dist;
    double a;
    double gamma;
    double sigma;
    double mass;
    
    bool pbc;
};

#endif	/* ARGUMENTS_H */