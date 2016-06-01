#include "MembraneTriangulation.h"

MembraneTriangulation::MembraneTriangulation() {}

MembraneTriangulation::MembraneTriangulation(const MembraneTriangulation& orig) {}

MembraneTriangulation::~MembraneTriangulation() {}

std::list<Triangle> MembraneTriangulation::triangulate()
{
    triangulate(1.0, 0.25, 5);
    return tris;
}

std::list<Triangle> MembraneTriangulation::triangulate(double L, double eps, int N)
{
//    N = 6;
//    double a = 2*(L + eps) / N;
//
//    put2Triangles(L + eps);
//    for (int i = 0; i < N; i++)
//    {
//        subdivide();
//    }

    double ta = (L + eps) / N ;
    double h = 0.5 * sqrt(3.0) * ta;
    createHex(ta, h);

    Vector3D x_shift(ta, 0, 0);
    Vector3D y_shift(0, 2 * h, 0);
    int Nx = (int) (L + eps) / (ta);
    int Ny = (int) (L + eps) / (2 * h);

    Vector3D shift;

    for (int ix = -Nx; ix < Nx; ix++)
    {
        for (int jy = -Ny; jy < Ny; jy++)
        {

            for (std::list<Triangle>::iterator i = hexagon.begin(); i != hexagon.end(); ++i)
            {
                shift = ix * x_shift + jy * y_shift;
            }

            for (std::list<Triangle>::iterator i = diamond1.begin(); i != diamond1.end(); ++i)
            {
                shift = ix * x_shift + jy * y_shift;
                tris.push_back(Triangle( shift + i->a, shift + i->b, shift + i->c ));
            }

            for (std::list<Triangle>::iterator i = diamond2.begin(); i != diamond2.end(); ++i)
            {
                shift = ix * x_shift + jy * y_shift;
                tris.push_back(Triangle( shift + i->a, shift + i->b, shift + i->c ));
            }
        }
    }


    return tris;

}

void MembraneTriangulation::createHex(double a, double h)
{
    double a2 = a / 2.0;

    hexagon.push_back(Triangle(Vector3D( 0, 0, 0), Vector3D(a, 0, 0), Vector3D( a2, +h, 0)));
    hexagon.push_back(Triangle(Vector3D( 0, 0, 0), Vector3D(a, 0, 0), Vector3D( a2, -h, 0)));

    hexagon.push_back(Triangle(Vector3D( 0, 0, 0), Vector3D(a2, -h, 0), Vector3D( -a2, -h, 0)));
    hexagon.push_back(Triangle(Vector3D( 0, 0, 0), Vector3D(a2, h, 0), Vector3D( -a2, h, 0)));

    hexagon.push_back(Triangle(Vector3D( 0, 0, 0), Vector3D(-a, 0, 0), Vector3D( -a2, h, 0)));
    hexagon.push_back(Triangle(Vector3D( 0, 0, 0), Vector3D(-a, 0, 0), Vector3D( -a2, -h, 0)));

    diamond1.push_back(Triangle(Vector3D( 0, 0, 0), Vector3D(a2, -h, 0), Vector3D( -a2, -h, 0)));
    diamond1.push_back(Triangle(Vector3D( 0, 0, 0), Vector3D(a2, h, 0), Vector3D( -a2, h, 0)));

    diamond2.push_back(Triangle(Vector3D( 0, 0, 0), Vector3D(a, 0, 0), Vector3D( a2, +h, 0)));
    diamond2.push_back(Triangle(Vector3D( 0, 0, 0), Vector3D(a, 0, 0), Vector3D( a2, -h, 0)));
}

void MembraneTriangulation::putTwoTriangles(double L0)
{
    tris.push_back(Triangle(Vector3D( L0, L0, 0), Vector3D( -L0, -L0, 0), Vector3D( L0, -L0, 0) ) );
    tris.push_back(Triangle(Vector3D( L0, L0, 0), Vector3D( -L0, -L0, 0), Vector3D( -L0, L0, 0) ) );
}

void MembraneTriangulation::subdivide()
{
    int counter = 0;
    std::list<Triangle> newTris;

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)  // go through all triangles
    {
        Vector3D mid = (i->a + i->b) * 0.5f; // point between points A and B
        newTris.push_back(Triangle(i->b, i->c, mid)); // remember new triangles
        newTris.push_back(Triangle(i->a, i->c, mid));
        counter++;
    }

    tris.swap(newTris); // use new set of triangles;
}