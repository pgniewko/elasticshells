#include "SurfaceForce.h"

SurfaceForce::SurfaceForce() {}

SurfaceForce::SurfaceForce(const SurfaceForce& orig) {}

SurfaceForce::~SurfaceForce() {}

double SurfaceForce::calcForces(Box& box, std::vector<Cell>& cells)
{
    if (box.pbc)
    {
        return 0.0;
    }

//    double ecw = box.E_box;
//    double rvertex;
    Vector3D wallYZ, wallXZ, wallXY;
    Vector3D vertXYZ;
    double sgnx, sgny, sgnz;
    double bsx = box.getX();
    double bsy = box.getY();
    double bsz = box.getZ();
    double fx, fy, fz;
    Vector3D forceX(0, 0, 0);
    Vector3D forceY(0, 0, 0);
    Vector3D forceZ(0, 0, 0);
    int numOfCells = cells.size();
    double totalForce = 0.0;
    Vector3D djk;

    double eb = box.getE();
    double rb_= 0.0;
    double nub = box.getNu();

    double e1;
    double r1;
    double nu1;

    for (int i = 0; i < numOfCells; i++)
    {
        //rvertex = cells[i].getVertexR();

        e1 = cells[i].getE();
        r1 = cells[i].getVertexR();
        nu1 = cells[i].getNu();

        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            vertXYZ = cells[i].getVertexXYZ(j);
            sgnx = SIGN(vertXYZ.x);
            wallYZ.x = sgnx * bsx;
            wallYZ.y = vertXYZ.y;
            wallYZ.z = vertXYZ.z;
            djk = vertXYZ - wallYZ;
            //forceX = HertzianRepulsion::calcForce(djk, rvertex, ecw);
            forceX = HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
            fx = forceX.length();
            sgny = SIGN(vertXYZ.y);
            wallXZ.x = vertXYZ.x;
            wallXZ.y = sgny * bsy;
            wallXZ.z = vertXYZ.z;
            djk = vertXYZ - wallXZ;
            //forceY = HertzianRepulsion::calcForce(djk, rvertex, ecw);
            forceY = HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
            fy = forceY.length();
            sgnz = SIGN(vertXYZ.z);
            wallXY.x = vertXYZ.x;
            wallXY.y = vertXYZ.y;
            wallXY.z = sgnz * bsz;
            djk = vertXYZ - wallXY;
            //forceZ = HertzianRepulsion::calcForce(djk, rvertex, ecw);
            forceZ = HertzianRepulsion::calcForce(djk, r1, rb_, e1, eb, nu1, nub);
            fz = forceZ.length();
            totalForce +=  (fx + fy + fz);
        }
    }

    return totalForce;
}
