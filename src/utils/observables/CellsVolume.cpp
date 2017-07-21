#include "CellsVolume.h"

CellsVolume::CellsVolume (const char* name, const char* format) : Observer(name, format) {}

CellsVolume::CellsVolume(const CellsVolume& orig) : Observer(orig) {}

CellsVolume::~CellsVolume() {}

void CellsVolume::set_params(const int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
}

double CellsVolume::observe(const Box& box, std::vector<Cell>& cells, const DomainList& dl)
{
    return calcCellsVolume(cells);
}

double CellsVolume::calcCellsVolume(std::vector<Cell>& cells)
{
    double cellsVolume = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        cellsVolume += cells[i].calcVolume(d_param);
    }

    if ( cells[0].getNumberVertices() == 1 )
    {
        double overlapvolume = 0.0;

        for (uint i = 0; i < cells.size(); i++)
        {
            for (uint j = i + 1; j < cells.size(); j++)
            {
                overlapvolume += sphereSphereIntersection(cells[i], cells[j]);
            }
        }

        cellsVolume -= overlapvolume;
    }

    return cellsVolume;
}

double CellsVolume::sphereSphereIntersection(const Cell& c1, const Cell& c2)
{

    double r1 = c1.getInitR();
    double r2 = c2.getInitR();
    Vector3D cm1 = c1.getCm();
    Vector3D cm2 = c2.getCm();

    double d = (cm1 - cm2).length();

    double V = 0.0;

    double rmin = std::min(r1, r2);

    if ( d > (r1 + r2) )
    {
        V = 0.0;
    }
    else if ( d <= std::abs(r1 - r2) )
    {
        V = 4.0 * constants::pi * rmin * rmin * rmin / 3.0;
    }
    else
    {
        V  = (r1 + r2 - d) * (r1 + r2 - d) * constants::pi;
        V *= (d * d + 2 * d * (r1 + r2) - 3 * (r1 * r1 + r2 * r2) + 6 * r1 * r2);
        V /= (12 * d);
    }

    return V;
}

DerivedRegister<CellsVolume> CellsVolume::reg("CellsVolume");