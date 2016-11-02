#include "RandomTriangulation.h"

RandomTriangulation::RandomTriangulation(int ns, int na, double tmin, double tmax, double rv) : n_steps(ns), n_anneals(na), T_min(tmin), T_max(tmax), r_vertex(rv)
{
}

RandomTriangulation::RandomTriangulation(const RandomTriangulation& orig) 
{
}

RandomTriangulation::~RandomTriangulation() 
{
}

std::list<Triangle> RandomTriangulation::triangulate(double r0)
{
    int n_ = 4.0 * (r0/r_vertex) * (r0/r_vertex);
    int nrow = 6;
    double scaled_sigma = 1.5 * r_vertex * (1.0 / r0);
    
    double* xyz = new double[3 * n_];
    int* ltri = new int[nrow * 2 * (n_ - 2)];
    
    generate_random_points(n_, xyz, n_steps, n_anneals, T_min, T_max, scaled_sigma);
    
    int v1_idx;
    int v2_idx;
    int v3_idx;
    
    for (int i = 0; i < nrow * 2 * (n_ - 2); i++)
    {
        v1_idx = ltri[6*i + 0];
        v2_idx = ltri[6*i + 1];
        v3_idx = ltri[6*i + 2];
        
        ri[6*i + 0] << " " << ltri[6*i + 1] << " "<< ltri[6*i + 2] << "||" <<  ltri[6*i + 3] << " " << ltri[6*i + 4] << " "<< ltri[6*i + 5]  << endl;
    }    
    
    
    delete[] xyz;
    delete[] ltri; 
    
}
