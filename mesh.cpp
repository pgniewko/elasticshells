/*
 * Author : Pawel Gniewek (UC Berkeley)
 * Email  : pawel.gniewek@berkeley.edu
 * License: BSD
 */

#include <iostream>    /* cout, cin */
#include <fstream>     /* ios, ofstream*/
#include <stdio.h>     /* printf, fgets */
#include <argp.h>      /* argp_parse */
#include <stdlib.h>    /* atoi,  strtod */
#include <math.h>      /* log, sqrt */
#include <string>
#include <limits>
#include <iomanip>

#include <stdio.h>
#include <string.h>

#include "Environment.h"
#include "src/Timer.h"

#include "utils/Logger.h"
#include "utils/LogManager.h"
#include "geometry/Triangle.h"
#include "geometry/algorithms/Triangulation.h"
#include "geometry/algorithms/SimpleTriangulation.h"
#include "geometry/algorithms/PlatonicTriangulatoin.h"
#include "geometry/algorithms/RandomTriangulation.h"
#include "src/Shell.h"

//Timer clocks;
//double simulation_time;

//bool isUnique(std::list<Vector3D>& vlist, Vector3D& v)
//{
//    for (std::list<Vector3D>::iterator i = vlist.begin(); i != vlist.end(); ++i)
//    {
//        if (i->x == v.x && i->y == v.y && i->z == v.z)
//        {
//            return false;
//        }
//    }
//
//    return true;
//}

bool isUnique(std::list<Vector3D>& vlist, Vector3D& v, double e)
{
    for (std::list<Vector3D>::iterator i = vlist.begin(); i != vlist.end(); ++i)
    {
        if ( fabs(i->x - v.x) < e && fabs(i->y - v.y) < e && fabs(i->z - v.z) < e )
        {
            return false;
        }
    }

    return true;
}

int constructVertices(std::list<Triangle>& tris, Vertex* vertices)
{
    std::list<Vector3D> vectors;
    double xtmp, ytmp, ztmp;

    int number_v = 0;

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        if ( isUnique(vectors, i->a, 0.001) )
        {
            vectors.push_back(i->a);
            xtmp = i->a.x;
            ytmp = i->a.y;
            ztmp = i->a.z;
            vertices[number_v] = Vertex(xtmp, ytmp, ztmp);
            vertices[number_v].set_id(number_v);
            number_v++;
        }

        if ( isUnique(vectors, i->b, 0.001) )
        {
            vectors.push_back(i->b);
            xtmp = i->b.x;
            ytmp = i->b.y;
            ztmp = i->b.z;
            vertices[number_v] = Vertex(xtmp, ytmp, ztmp);
            vertices[number_v].set_id(number_v);
            number_v++;
        }

        if ( isUnique(vectors, i->c, 0.001) )
        {
            vectors.push_back(i->c);
            xtmp = i->c.x;
            ytmp = i->c.y;
            ztmp = i->c.z;
            vertices[number_v] = Vertex(xtmp, ytmp, ztmp);
            vertices[number_v].set_id(number_v);
            number_v++;
        }
    }

    return number_v;
}

int getVertex(Vector3D& v, Vertex* vertices, int number_v, double e)
{
    for (int i = 0; i < number_v; i++)
    {
        if ( fabs(vertices[i].r_c.x - v.x)  < e && fabs(vertices[i].r_c.y - v.y) < e && fabs(vertices[i].r_c.z - v.z) < e)
        {
            return vertices[i].get_id();
        }
    }

    return -1;
}


int  constructVTriangles(std::list<Triangle>& tris, Vertex* vertices, VertexTriangle* triangles, int number_v)
{
    int number_t = 0;

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        int va = getVertex(i->a, vertices, number_v, 0.001);
        int vb = getVertex(i->b, vertices, number_v, 0.001);
        int vc = getVertex(i->c, vertices, number_v, 0.001);
        VertexTriangle vrxt(va, vb, vc);
        triangles[number_t] = VertexTriangle(vrxt);
        triangles[number_t].setId(number_t);
        number_t++;
    }

    return number_t;
}




double min_theta(Triangle tri)
{
    Vector3D ab = tri.a - tri.b;
    Vector3D ac = tri.a - tri.c;
    Vector3D bc = tri.b - tri.c;

    double a = ab.length();
    double b = ac.length();
    double c = bc.length();

    double cosA = std::cos(a * a + b * b - c * c) / (2 * a * b);
    double cosB = std::cos(a * a + c * c - b * b) / (2 * a * c);
    double cosC = std::cos(b * b + c * c - a * a) / (2 * b * c);

    double thetaA = std::acos(cosA);
    double thetaB = std::acos(cosB);
    double thetaC = std::acos(cosC);

    double mintheta = std::numeric_limits<double>::max();
    mintheta = std::min(mintheta, thetaA);
    mintheta = std::min(mintheta, thetaB);
    mintheta = std::min(mintheta, thetaC);
    return mintheta;
}

double q(std::list<Triangle>& tris)
{
    double min_q = std::numeric_limits<double>::max();

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        Vector3D ab = i->a - i->b;
        Vector3D ac = i->a - i->c;
        Vector3D bc = i->b - i->c;

        double a = ab.length();
        double b = ac.length();
        double c = bc.length();
        min_q = std::min( min_q, (b + c - a) * (c + a - b) * (a + b - c) / (a * b * c));
    }


    return min_q;
}

double q_A(Vertex* vertices, VertexTriangle* triangles, int number_v)
{
    int idx;
    double min_A = std::numeric_limits<double>::max();
    double max_A = std::numeric_limits<double>::min();

    double tmp;

    for (int i = 0; i < number_v; i++)
    {
        tmp = 0.0;

        for (int j = 0; j < vertices[i].facets_number; j++)
        {
            idx = vertices[i].bondedTris[j];
            tmp += triangles[idx].area(vertices) / 3.0;
        }

        min_A = std::min(min_A, tmp);
        max_A = std::max(max_A, tmp);
    }

    return min_A / max_A;
}

double q_theta(std::list<Triangle>& tris)
{
    double min_angle = std::numeric_limits<double>::max();

    double tmp;

    for (std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
        tmp = i->min_angle();
        min_angle = std::min(min_angle, tmp);
    }

    return min_angle / (M_PI / 3.0);
}

void constructTopology(Vertex* vertices, VertexTriangle* triangles, int number_v, int number_t)
{
    for (int i = 0; i < number_t; i++)
    {
        int aid = triangles[i].ia;
        int bid = triangles[i].ib;
        int cid = triangles[i].ic;
        Vector3D ab = vertices[aid].r_c - vertices[bid].r_c;
        Vector3D ac = vertices[aid].r_c - vertices[cid].r_c;
        Vector3D bc = vertices[bid].r_c - vertices[cid].r_c;
        double abl = ab.length();
        double acl = ac.length();
        double bcl = bc.length();
        int tid = triangles[i].my_id;
        vertices[aid].addNeighbor(bid, abl);
        vertices[aid].addNeighbor(cid, acl);
        vertices[bid].addNeighbor(aid, abl);
        vertices[bid].addNeighbor(cid, bcl);
        vertices[cid].addNeighbor(aid, acl);
        vertices[cid].addNeighbor(bid, bcl);
        vertices[aid].addTriangle(tid);
        vertices[bid].addTriangle(tid);
        vertices[cid].addTriangle(tid);
    }
}

int main(int argc, char** argv)
{

    Vertex* vertices = new Vertex[170000];
    VertexTriangle* triangles = new VertexTriangle[330000];

    int number_v = 0;
    int number_t = 0;

    //    int depths[8] = {1, 2, 3, 4, 5, 6, 7, 8};
//    std::cout << "Level & Points & Triangles & $q_A$ & $q_\\theta$ & $q$ \\\\  [0.5ex]  \\hline \\hline" << std::endl;
//    for (int i = 0; i < 8; i++)
//    {
//        SimpleTriangulation mesh(i + 1);
//        //PlatonicTriangulatoin mesh(i+1, 3);
//        std::list<Triangle> tris = mesh.triangulate();
//        number_v = constructVertices(tris, vertices);
//        number_t = constructVTriangles(tris, vertices, triangles, number_v);
//        constructTopology(vertices, triangles, number_v, number_t);
//
//        std::cout << depths[i] << " & " << number_v << " & " << tris.size() << " & " << std::setprecision(5) << q_A(vertices, triangles, number_v) << " & " << std::setprecision(5) << q_theta(tris) << " & " << std::setprecision(5) << q(tris) << " \\\\ [1ex] \\hline" << std::endl;
//    }


    double rvs[] = {0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1};
    std::cout << "Rvs & Points & Triangles & $q_A$ & $q_\\theta$ & $q$ \\\\  [0.5ex]  \\hline \\hline" << std::endl;

    double qa, qt, qr;

    for (int i = 0; i < 9; i++)
    {
        qa = 0.0;
        qt = 0.0;
        qr = 0.0;

        for (int j = 0; j < 5; j++)
        {
            RandomTriangulation rnd(100, 200, 0.1, 1000.0, rvs[i]);
            std::list<Triangle> tris = rnd.triangulate(2.5);

            number_v = constructVertices(tris, vertices);
            number_t = constructVTriangles(tris, vertices, triangles, number_v);
            constructTopology(vertices, triangles, number_v, number_t);
            qa += q_A(vertices, triangles, number_v);
            qt += q_theta(tris);
            qr  += q(tris);

        }

        qa /= 5;
        qt /= 5;
        qr /= 5;
        std::cout << rvs[i] << " & " << number_v << " & " << number_t << " & " << std::setprecision(5) << qa << " & " << std::setprecision(5) << qt << " & " << std::setprecision(5) << qr << " \\\\ [1ex] \\hline" << std::endl;
    }

    double temps[] = {0.1, 1.0, 10, 50, 100, 250, 500, 750};
    std::cout << "Temp & Points & Triangles & $q_A$ & $q_\\theta$ & $q$ \\\\  [0.5ex]  \\hline \\hline" << std::endl;

    for (int i = 0; i < 8; i++)
    {
        qa = 0.0;
        qt = 0.0;
        qr = 0.0;

        for (int j = 0; j < 5; j++)
        {
            RandomTriangulation rnd(100, 200, temps[i], 1000.0, 0.1);
            std::list<Triangle> tris = rnd.triangulate(2.5);

            number_v = constructVertices(tris, vertices);
            number_t = constructVTriangles(tris, vertices, triangles, number_v);
            constructTopology(vertices, triangles, number_v, number_t);
            qa += q_A(vertices, triangles, number_v);
            qt += q_theta(tris);
            qr += q(tris);
        }

        qa /= 5;
        qt /= 5;
        qr /= 5;
        std::cout << temps[i] << " & " << number_v << " & " << number_t << " & " << std::setprecision(5) << qa << " & " << std::setprecision(5) << qt << " & " << std::setprecision(5) << qr << " \\\\ [1ex] \\hline" << std::endl;

    }




    delete[] vertices;
    delete[] triangles;

    return 0;
}

