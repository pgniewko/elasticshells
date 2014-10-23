#include "SurfacePressure.h"

SurfacePressure::SurfacePressure() {}

SurfacePressure::SurfacePressure(const SurfacePressure& orig) {}

SurfacePressure::~SurfacePressure() {}

double SurfacePressure::calcPressure(Box& box, vector<Cell>& cells)
{
    double ecs = box.ecw;
    double rcb;
    Vector3D wallYZ, wallXZ, wallXY;
    Vector3D vertXYZ;
    
    double sgnx, sgny, sgnz;
    double bsx = box.getX();
    double bsy = box.getY();
    double bsz = box.getZ();

    Vector3D forceX(0, 0, 0);
    Vector3D forceY(0, 0, 0);
    Vector3D forceZ(0, 0, 0);
    
    double area = box.getArea();
    int numOfCells = cells.size();
    double totalForce = 0.0;
    Vector3D djk;
    
    for (int i = 0; i < numOfCells; i++)
    {
        rcb = cells[i].getRcb();
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            vertXYZ = cells[i].getVertexXYZ(j);
            
            if (vertXYZ.x != 0 )
            {
                sgnx = vertXYZ.x / fabs(vertXYZ.x);
                wallYZ.x = sgnx * bsx;
                wallYZ.y = vertXYZ.y;
                wallYZ.z = vertXYZ.z;
                djk = vertXYZ - wallYZ;
                forceX = BoxCellRepulsion::calcForce(djk, ecs, rcb);
                totalForce += forceX.length();
            }

            if (vertXYZ.y != 0 )
            {
                sgny = vertXYZ.y / fabs(vertXYZ.y);
                wallXZ.x = vertXYZ.x;
                wallXZ.y = sgny * bsy;
                wallXZ.z = vertXYZ.z;
                djk = vertXYZ - wallXZ;
                forceY = BoxCellRepulsion::calcForce(djk, ecs, rcb);
                totalForce += forceY.length();
            }

            if (vertXYZ.z != 0 )
            {
                sgnz = vertXYZ.z / fabs(vertXYZ.z);
                wallXY.x = vertXYZ.x;
                wallXY.y = vertXYZ.y;
                wallXY.z = sgnz * bsz;
                djk = vertXYZ - wallXY;
                forceZ = BoxCellRepulsion::calcForce(djk, ecs, rcb);
                totalForce += forceZ.length();
            }
        }
        
        //forcemagn += cells[i].calcBoxForces(box);
    }
    //double pressure = totalForce / area;
    return totalForce / area;
}
