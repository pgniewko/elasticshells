#include "Aspherity.h"

Aspherity::Aspherity(const char* name, const char* format) : Observer(name, format) {}

Aspherity::Aspherity(const Aspherity& orig) : Observer(orig) {}

Aspherity::~Aspherity() {}

void Aspherity::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double Aspherity::observe(const Box& box, const std::vector<Shell>& shells)
{
    double av_radius = 0.0;
    double sq_av_radius = 0.0;
    double sq_sum = 0.0;
    double sq_total = 0.0;

    for (uint i = 0; i < shells.size(); i++)
    {
        //shells[i].calc_cm();
        Vector3D cell_cm = shells[i].center_of_mass;

        for (int j = 0; j < shells[i].getNumberVertices(); j++)
        {
            av_radius += (shells[i].vertices[j].r_c - cell_cm).length();
        }

        av_radius /= shells[i].getNumberVertices();
        sq_av_radius = av_radius * av_radius;
        double res;

        for (int j = 0; j < shells[i].getNumberVertices(); j++)
        {
            res = (shells[i].vertices[j].r_c - cell_cm).length() - av_radius;
            sq_sum += res * res;
        }

        sq_sum /= sq_av_radius;
        sq_sum /= shells[i].getNumberVertices();
        sq_total += sq_sum;
    }

    sq_total /= shells.size();
    return sq_total;
}

DerivedRegister<Aspherity> Aspherity::reg("Aspherity");