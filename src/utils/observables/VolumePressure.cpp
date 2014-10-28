#include "VolumePressure.h"

VolumePressure::VolumePressure() {}

VolumePressure::VolumePressure(const VolumePressure& orig) {}

VolumePressure::~VolumePressure() {}

double VolumePressure::calcPressure(Box& box, std::vector<Cell>& cells, int nbhandler)
{
    int numberofCells = cells.size();
    double V = box.getVolume();
    double sab[] = {0, 0, 0, 0, 0, 0};
    
    if (!box.pbc)
    {
        // RESET FORCES
        for (int i = 0 ; i < numberofCells; i++)
        {
            cells[i].voidForces();
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
    }
    
    Vector3D vertXYZ;
    Vector3D vertForce;
    
    for (int i = 0; i < numberofCells; i++)
    {
        for (int j = 0; j < cells[i].numberOfVerts(); j++)
        {
            vertXYZ = cells[i].getVertexXYZ(j);
            vertForce = cells[i].getVertexForce(j);
            
            sab[0] += -vertXYZ.x  * vertForce.x;
            sab[1] += -vertXYZ.y  * vertForce.y;
            sab[2] += -vertXYZ.z  * vertForce.z;
            sab[3] += -vertXYZ.x  * vertForce.y;
            sab[4] += -vertXYZ.x  * vertForce.z;
            sab[5] += -vertXYZ.y  * vertForce.z;
        }
    
    }
    
    if (!box.pbc)
    {
        // RESET FORCES
        for (int i = 0 ; i < numberofCells; i++)
        {
            cells[i].voidForces();
        }
    }
    
    double pressure = -(sab[0] + sab[1] + sab[2]);
    pressure /= 3.0;
    pressure /= V;
    return pressure;
}