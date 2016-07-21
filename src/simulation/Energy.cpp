#include "Energy.h"

utils::Logger Energy::energy_logs("energy");
ulong ENERGY_EVALUATION_COUNTER(0);

Energy::Energy() {}

Energy::Energy(const Energy& orig) {}

Energy::~Energy() {}

double Energy::calcTotalEnergy(const std::vector<Cell>& cells, const Box& box)
{
    
    double totalEnergy  = 0.0;
    
    totalEnergy += calcMembraneEnergy(cells, box);
    totalEnergy += calcOsmoticEnergy(cells, box);
    totalEnergy += calcContactEnergy(cells, box);
    
    return totalEnergy;

}

double Energy::calcMembraneEnergy(const std::vector<Cell>& cells, const Box& box)
{
    double membraneEnergy = 0.0;
    membraneEnergy += calcStretchEnergy(cells, box);
    membraneEnergy += calcBendingEnergy(cells);
    
    return membraneEnergy;
}

double Energy::calcStretchEnergy(const std::vector<Cell>& cells, const Box& box, char* model_t)
{
    double stretchEnergy = 0.0;
    
    if (STRCMP (model_t, "fem"))
    {
        for (uint i = 0; i < cells.size(); i++)
        {
            for (int j = 0; j < cells[i].number_t; j++)
            {    
                stretchEnergy += cells[i].triangles[j].calcFemEnergy(cells[i].vertices);
            }
        }
    }
    else
    {
        // EXIT
    }
    return stretchEnergy;
}

double Energy::calcBendingEnergy(const std::vector<Cell>& cells)
{
    double bendingEnergy = 0.0;
    
    for (uint i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].number_s; j++)
        {
            bendingEnergy += cells[i].bhinges[j].calcBendingEnergy(cells[i].vertices);
        }
    }
    
    return bendingEnergy;
}

double Energy::calcOsmoticEnergy(const std::vector<Cell>& cells)
{
    double osmoticEnergy = 0.0;
    
    double P, V, Ve, nRT;
    for (uint i = 0; i < cells.size(); i++)
    {
        if ( OsmoticForce::getFlag() )
        {
            nRT = cells[i].nRT;
            Ve = cells[i].V0 * OsmoticForce::getEpsilon();
            V = cells[i].calcVolume();
            osmoticEnergy += -nRT*std::log(V - Ve);
        }
        else
        {
            P = cells[i].getTurgor();
            V = cells[i].calcVolume();
            osmoticEnergy += -P*V;
        }
    }
    
    return osmoticEnergy;
}

