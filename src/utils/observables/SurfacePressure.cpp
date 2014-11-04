#include "SurfacePressure.h"

SurfacePressure::SurfacePressure() {}

SurfacePressure::SurfacePressure(const SurfacePressure& orig) {}

SurfacePressure::~SurfacePressure() {}

double SurfacePressure::calcPressure(Box& box, std::vector<Cell>& cells)
{
    if (box.pbc)
        return 0.0;
    
    double maxRbc = 0.0;
    double ecs = box.ecw;
    double rcb;
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
    
    for (int i = 0; i < numOfCells; i++)
    {
        rcb = cells[i].getRbc();
        maxRbc = std::max(maxRbc, rcb);
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
        {
            vertXYZ = cells[i].getVertexXYZ(j);
           
            sgnx = SIGN(vertXYZ.x);
            wallYZ.x = sgnx * bsx;
            wallYZ.y = vertXYZ.y;
            wallYZ.z = vertXYZ.z;
            djk = vertXYZ - wallYZ;
            forceX = HertzianRepulsion::calcForce(djk, rcb, ecs);
            fx = forceX.length();
            
            sgny = SIGN(vertXYZ.y);
            wallXZ.x = vertXYZ.x;
            wallXZ.y = sgny * bsy;
            wallXZ.z = vertXYZ.z;
            djk = vertXYZ - wallXZ;
            forceY = HertzianRepulsion::calcForce(djk, rcb, ecs);
            fy = forceY.length();
            
            sgnz = SIGN(vertXYZ.z);
            wallXY.x = vertXYZ.x;
            wallXY.y = vertXYZ.y;
            wallXY.z = sgnz * bsz;
            djk = vertXYZ - wallXY;
            forceZ = HertzianRepulsion::calcForce(djk, rcb, ecs);
            fz = forceZ.length();
            totalForce +=  sqrt( fx*fx + fy*fy + fz*fz );
        }
    }
    
    //double area = box.getArea(maxRbc);
    double area = box.getArea();
    return totalForce / area;
}
