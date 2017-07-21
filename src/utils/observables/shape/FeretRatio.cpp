#include "FeretRatio.h"

FeretRatio::FeretRatio(const char* name, const char* format) : Observer(name, format) {}

FeretRatio::FeretRatio(const FeretRatio& orig) : Observer(orig) {}

FeretRatio::~FeretRatio() {}

void FeretRatio::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double FeretRatio::observe(const Box& box, std::vector<Cell>& cells, const DomainList& dl)
{
    double max_d = std::numeric_limits<double>::min();
    double D = 0.0, V;
    double FR = 0.0;


    for (uint i = 0; i < cells.size(); i++)
    {
        max_d = std::numeric_limits<double>::min();

        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            for (int k = j + 1; k < cells[i].getNumberVertices(); k++)
            {
                // MAXIMUM FERET'S DIAMETER a.k.a. Maximum Caliber
                max_d = std::max(max_d, (cells[i].vertices[j].r_c - cells[i].vertices[k].r_c).length() );
            }
        }

        V = cells[i].calcVolume();
        V = 3.0 * V / (4.0 * constants::pi);
        D = 2.0 * std::pow( V, 1.0 / 3.0 ); // HEYWOOD'S DIAMETER

        if (D > 0)
        {
            FR += (max_d / D) - 1.0;
        }
        else
        {
            throw DataException("Division by zero");
        }
    }

    FR /= cells.size();

    return FR;
}

DerivedRegister<FeretRatio> FeretRatio::reg("FeretRatio");