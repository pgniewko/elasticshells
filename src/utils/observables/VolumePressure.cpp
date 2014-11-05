#include "VolumePressure.h"

VolumePressure::VolumePressure() {}

VolumePressure::VolumePressure(const VolumePressure& orig) {}

VolumePressure::~VolumePressure() {}

double VolumePressure::calcPressure(Box& box, std::vector<Cell>& cells, int nbhandler)
{
    
    int numberofCells = cells.size();
    double V;// = box.getVolume();
    double sab[] = {0, 0, 0, 0, 0, 0};
    double maxRbc = 0.0;
    
    // RESET FORCES
    for (int i = 0 ; i < numberofCells; i++)
    {
        cells[i].voidForces();
        maxRbc = std::max(maxRbc, cells[i].getRbc());
    }

    // CALCULATE INTRA-CELLULAR FORCES
    for (int i = 0 ; i < numberofCells; i++)
    {
        cells[i].calcBondedForces();
    }

    // CALCULATE INTER-CELLULAR FORCES
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < numberofCells; j++)
        {
            if (nbhandler == 0)
            {
                cells[i].calcNbForcesON2(cells[j], box);
            }
            else if (nbhandler == 1)
            {
                cells[i].calcNbForcesVL(cells[j], box);
            }   
            else if (nbhandler == 2)
            {
                cells[i].calcNbForcesVL(cells[j], box);
            }
            else 
            {
                cells[i].calcNbForcesON2(cells[j], box);
            }
        }
    }
    
    //Vector3D vertXYZ;
    Vector3D VertXYZimg;
    Vector3D vertForce;
    
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
        {
            //vertXYZ = cells[i].getVertexXYZ(j);
            VertXYZimg = getImage(box, cells[i].getVertexXYZ(j) );
            vertForce = cells[i].getVertexForce(j);
            
            sab[0] += -VertXYZimg.x  * vertForce.x;
            sab[1] += -VertXYZimg.y  * vertForce.y;
            sab[2] += -VertXYZimg.z  * vertForce.z;
            sab[3] += -VertXYZimg.x  * vertForce.y;
            sab[4] += -VertXYZimg.x  * vertForce.z;
            sab[5] += -VertXYZimg.y  * vertForce.z;
        }
    
    }
    

    // RESET FORCES
    for (int i = 0 ; i < numberofCells; i++)
    {
        cells[i].voidForces();
    }
    
    //V = box.getVolume(maxRbc);
    V = box.getVolume();
    double pressure = -(sab[0] + sab[1] + sab[2]);
    pressure /= 3.0;
    pressure /= V;
    return pressure;
}

Vector3D VolumePressure::getImage(Box& box, const Vector3D& vec)
{
    Vector3D image(vec);
    if (!box.pbc)
        return image;
    
    
    double bsx = 2 * box.getX();
    double bsy = 2 * box.getY();
    double bsz = 2 * box.getZ();
    
    image.x -= bsx * rint( image.x / bsx );
    image.y -= bsy * rint( image.y / bsy );
    image.z -= bsz * rint( image.z / bsz );
    return image;
}