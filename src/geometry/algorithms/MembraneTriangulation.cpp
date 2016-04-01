#include "MembraneTriangulation.h"

MembraneTriangulation::MembraneTriangulation() {
}

MembraneTriangulation::MembraneTriangulation(const MembraneTriangulation& orig) {}

MembraneTriangulation::~MembraneTriangulation() {}

std::list<Triangle> MembraneTriangulation::triangulate()
{
    triangulate(1.0, 0.25, 5);
    return tris;
}

std::list<Triangle> MembraneTriangulation::triangulate(double L, double eps, int N)
{
    double a = 2*(L + eps) / N;
    double h = 0.5 * sqrt(3.0) * a;
    createHex(L, eps, a);
    Vector3D x_shift(a,0,0);
    Vector3D y_shift(0,2*h,0);
    
    int Nx = (int) 2*(L + eps) / (a);
    int Ny = (int) 2*(L + eps) / (2*h);
    
    std::cout << "Nx="<<Nx<< std::endl;
    std::cout << "Ny="<<Ny<< std::endl;
    
    for (int ix = -Nx; ix <= Nx; ix++)
    {
        for (int jy = -Ny; jy <= Ny; jy++)
        {
            if (ix!=0 || jy!=0)
            {
                std::cout << "ix="<<ix << " iy=" << jy<< std::endl;
                for (std::list<Triangle>::iterator i = hexagon.begin(); i != hexagon.end(); ++i)
                {
                    Vector3D shift = ix*x_shift+ + jy*y_shift;
                    tris.push_back(Triangle( shift+i->a, shift+i->b, shift+i->c ));
                }
            }
        }    
    }
    
    tris.push_back( Triangle( Vector3D( 0,0,0 ), Vector3D( 0,0,0 ), Vector3D (0,0,0 ) ) ) ;
     
    return tris;

}

void MembraneTriangulation::createHex(double L, double eps, double a)
{
    double s = 0.0;//-(L+eps);
    double a2 = a / 2.0;
    double h = 0.5 * sqrt(3.0) * a;
    
    hexagon.push_back(Triangle(Vector3D( 0, s, 0), Vector3D(s+a, s, 0), Vector3D( s+a2, s+h, 0)));
    hexagon.push_back(Triangle(Vector3D( s, s, 0), Vector3D(s+a, s, 0), Vector3D( s+a2, s-h, 0)));

    hexagon.push_back(Triangle(Vector3D( s, s, 0), Vector3D(s+a2, s-h, 0), Vector3D( s-a2, s-h, 0)));
    hexagon.push_back(Triangle(Vector3D( s, s, 0), Vector3D(s+a2, s+h, 0), Vector3D( s-a2, s+h, 0)));
 
    hexagon.push_back(Triangle(Vector3D( s, s, 0), Vector3D(s-a, s, 0), Vector3D( s-a2, s+h, 0)));
    hexagon.push_back(Triangle(Vector3D( s, s, 0), Vector3D(s-a, s, 0), Vector3D( s-a2, s-h, 0)));
    
    for (std::list<Triangle>::iterator i = hexagon.begin(); i != hexagon.end(); ++i)
    {
        tris.push_back(Triangle(i->a, i->b, i->c ));
    }
}