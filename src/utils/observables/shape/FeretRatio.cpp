#include "FeretRatio.h"

FeretRatio::FeretRatio(const char* name, const char* format) : Observer(name, format) {}

FeretRatio::FeretRatio(const FeretRatio& orig) : Observer(orig) {}

FeretRatio::~FeretRatio() {}

void FeretRatio::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double FeretRatio::observe(const Box& box, const std::vector<Shell>& shells)
{
    double max_d = 0.0;
    double D = 0.0, V;
    double FR = 0.0;


    for (uint i = 0; i < shells.size(); i++)
    {
        max_d = 0.0;

        for (int j = 0; j < shells[i].get_number_vertices(); j++)
        {
            for (int k = j + 1; k < shells[i].get_number_vertices(); k++)
            {
                // MAXIMUM FERET'S DIAMETER a.k.a. Maximum Caliber
                max_d = std::max(max_d, (shells[i].vertices[j].r_c - shells[i].vertices[k].r_c).length() );
            }
        }

        V = shells[i].calc_volume();
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

    FR /= shells.size();

    return FR;
}

DerivedRegister<FeretRatio> FeretRatio::reg("FeretRatio");