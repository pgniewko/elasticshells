#include "Elongation.h"

Elongation::Elongation(const char* name, const char* format) : Observer(name, format) {}

Elongation::Elongation(const Elongation& orig) : Observer(orig) {}

Elongation::~Elongation() {}

void Elongation::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double Elongation::observe(const Box& box, std::vector<Cell>& cells, const DomainList& dl)
{
    double min_r = std::numeric_limits<double>::max();
    double max_r = std::numeric_limits<double>::min();
    
    double el = 0.0;
    double l;

    for (uint i = 0; i < cells.size(); i++)
    {
        min_r = std::numeric_limits<double>::max();
        max_r = std::numeric_limits<double>::min();
        
        cells[i].calcCM();
        Vector3D cell_cm = cells[i].cm_m;

        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            l = (cells[i].vertices[j].r_c - cell_cm).length();
            min_r = std::min( min_r, l );
            max_r = std::max( max_r, l );
        }

        el += ( (max_r - min_r)/ min_r );
    }
    
    el /= cells.size();

}

DerivedRegister<Elongation> Elongation::reg("Elongation");

