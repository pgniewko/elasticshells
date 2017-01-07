#include "RandomTriangulation.h"

RandomTriangulation::RandomTriangulation(int ns, int na, double tmin, double tmax, double rv) : n_steps(ns), n_anneals(na), T_min(tmin), T_max(tmax), r_vertex(rv)
{
}

RandomTriangulation::RandomTriangulation(const RandomTriangulation& orig) : 
n_steps(orig.n_steps), n_anneals(orig.n_anneals), T_min(orig.T_min), T_max(orig.T_max), r_vertex(orig.r_vertex)
{}

RandomTriangulation::~RandomTriangulation() 
{}

std::list<Triangle> RandomTriangulation::triangulate()
{
    return triangulate(1.0);
}
std::list<Triangle> RandomTriangulation::triangulate(double r0)
{
    int n_ = 4.0 * (r0 / r_vertex) * (r0 / r_vertex);
    int nrow = 6;
    
    //std::cout << "RANDOM TRIANGUILATION; n_ =" << n_ << " r0=" << r0 << std::endl;
    
    double* xyz = new double[3 * n_];
    int* ltri = new int[nrow * 2 * (n_ - 2)];
    
    generate_random_points(n_, xyz, n_steps, n_anneals, T_min, T_max, r_vertex);
    
    traingulate_points(n_, xyz, ltri);
    
    int v1_idx;
    int v2_idx;
    int v3_idx;
    
    int tidx;
    
    for (int i = 0; i < 2 * (n_ - 2); i++)
    {

        tidx = i;
        v1_idx = ltri[nrow * i + 0]-1;
        v2_idx = ltri[nrow * i + 1]-1;
        v3_idx = ltri[nrow * i + 2]-1;
        
        Vector3D va = Vector3D(xyz[3*v1_idx + 0], xyz[3*v1_idx + 1], xyz[3*v1_idx + 2]);
        Vector3D vb = Vector3D(xyz[3*v2_idx + 0], xyz[3*v2_idx + 1], xyz[3*v2_idx + 2]);
        Vector3D vc = Vector3D(xyz[3*v3_idx + 0], xyz[3*v3_idx + 1], xyz[3*v3_idx + 2]);
        
        tris.push_back( Triangle(va, vb, vc) );
    }
    
    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        i->a.set_length(r0);
        i->b.set_length(r0);
        i->c.set_length(r0);
    }
    
    delete[] xyz;
    delete[] ltri; 
    
    return tris;
}
