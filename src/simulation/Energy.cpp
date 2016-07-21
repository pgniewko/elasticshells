#include "Energy.h"

utils::Logger Energy::energy_logs("energy");
ulong Energy::ENERGY_EVALUATION_COUNTER(0);

Energy::Energy() {}

Energy::Energy(const Energy& orig) {}

Energy::~Energy() {}

double Energy::calcTotalEnergy(const std::vector<Cell>& cells, const Box& box, const DomainList& domains, char* model_t)
{
    double totalEnergy  = 0.0;
    
    totalEnergy += calcMembraneEnergy(cells, box, model_t);
    totalEnergy += calcOsmoticEnergy(cells);
    totalEnergy += calcContactEnergy(cells, box, domains);
    
    ENERGY_EVALUATION_COUNTER++;
    
    return totalEnergy;
}

double Energy::calcMembraneEnergy(const std::vector<Cell>& cells, const Box& box, char* model_t)
{
    double membraneEnergy = 0.0;
    membraneEnergy += calcStretchEnergy(cells, box, model_t);
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
        if (!cells[i].no_bending)
        {
            for (int j = 0; j < cells[i].number_s; j++)
            {
                bendingEnergy += cells[i].bhinges[j].calcBendingEnergy(cells[i].vertices);
            }
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

double Energy::calcContactEnergy(const std::vector<Cell>& cells, const Box& box, const DomainList& domains)
{
    double totalContactEnergy = 0.0;
    totalContactEnergy += calcCellBoxEnergy(cells, box);
    totalContactEnergy += calcCellCellEnergy(cells, box, domains);
    
    return totalContactEnergy;
}

double Energy::calcCellBoxEnergy(const std::vector<Cell>& cells, const Box& box)
{
    double cell_box_energy = 0.0;
    
    
    Vector3D wallYZ, wallXZ, wallXY;
    Vector3D dij;
    double sgnx, sgny, sgnz;
    double bsx = box.getX();
    double bsy = box.getY();
    double bsz = box.getZ();
    double eb  = box.getE();
    double nub = box.getNu();
    double rb_ = 0.0;
    double r1;
    double e1;
    double nu1;
    
    for (uint i = 0; i < cells.size(); i++)
    {
        r1 = cells[i].params.vertex_r;
        e1 = cells[i].params.ecc;
        nu1 = cells[i].params.nu;
        
        for (int j = 0; j < cells[i].number_v; j++)
        {
            sgnx = SIGN(cells[i].vertices[j].r_c.x);
            wallYZ.x = sgnx * bsx;
            wallYZ.y = cells[i].vertices[j].r_c.y;
            wallYZ.z = cells[i].vertices[j].r_c.z;
            dij = cells[i].vertices[j].r_c - wallYZ;
            cell_box_energy += HertzianRepulsion::calcEnergy(dij, r1, rb_, e1, eb, nu1, nub);
            
            sgny = SIGN(cells[i].vertices[j].r_c.y);
            wallXZ.x = cells[i].vertices[j].r_c.x;
            wallXZ.y = sgny * bsy;
            wallXZ.z = cells[i].vertices[j].r_c.z;
            dij = cells[i].vertices[j].r_c - wallXZ;
            cell_box_energy += HertzianRepulsion::calcEnergy(dij, r1, rb_, e1, eb, nu1, nub);
            
            sgnz = SIGN(cells[i].vertices[j].r_c.z);
            wallXY.x = cells[i].vertices[j].r_c.x;
            wallXY.y = cells[i].vertices[j].r_c.y;
            wallXY.z = sgnz * bsz;
            dij = cells[i].vertices[j].r_c - wallXY;
            cell_box_energy += HertzianRepulsion::calcEnergy(dij, r1, rb_, e1, eb, nu1, nub);
        }   
    }
    
    return cell_box_energy;
}

double Energy::calcCellCellEnergy(const std::vector<Cell>& cells, const Box& box, const DomainList& domains)
{
    return domains.calcNbEnergy(cells, box);
}