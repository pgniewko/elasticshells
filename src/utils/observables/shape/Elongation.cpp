#include "Elongation.h"

Elongation::Elongation(const char* name, const char* format) : Observer(name, format) {}

Elongation::Elongation(const Elongation& orig) : Observer(orig) {}

Elongation::~Elongation() {}

double Elongation::observe(const Box& box, const std::vector<Shell>& shells)
{
    double min_r = 1000000.0;
    double max_r = 0.0;

    double el = 0.0;
    double l;

    for (uint i = 0; i < shells.size(); i++)
    {
        min_r = 1000000.0;
        max_r = 0.0;

        Vector3D cell_cm = shells[i].center_of_mass;

        for (int j = 0; j < shells[i].get_number_vertices(); j++)
        {
            l = (shells[i].vertices[j].r_c - cell_cm).length();
            min_r = std::min( min_r, l );
            max_r = std::max( max_r, l );
        }

        if (min_r > 0)
        {
            el += ( max_r / min_r  - 1.0 );
        }
        else
        {
            throw DataException("Division by zero");
        }
    }

    el /= shells.size();

    return el;

}

DerivedRegister<Elongation> Elongation::reg("Elongation");

