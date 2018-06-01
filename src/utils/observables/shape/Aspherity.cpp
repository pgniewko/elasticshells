#include "Aspherity.h"

Aspherity::Aspherity(const char* name, const char* format) : Observer(name, format) {}

Aspherity::Aspherity(const Aspherity& orig) : Observer(orig) {}

Aspherity::~Aspherity() {}

void Aspherity::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double Aspherity::observe(const Box& box, std::vector<Shell>& cells, const DomainList& dl)
{
    double av_radius = 0.0;
    double sq_av_radius = 0.0;
    double sq_sum = 0.0;
    double sq_total = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        cells[i].calcCM();
        Vector3D cell_cm = cells[i].cm_m;

        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            av_radius += (cells[i].vertices[j].r_c - cell_cm).length();
        }

        av_radius /= cells[i].getNumberVertices();
        sq_av_radius = av_radius * av_radius;
        double res;

        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            res = (cells[i].vertices[j].r_c - cell_cm).length() - av_radius;
            sq_sum += res * res;
        }

        sq_sum /= sq_av_radius;
        sq_sum /= cells[i].getNumberVertices();
        sq_total += sq_sum;
    }

    sq_total /= cells.size();
    return sq_total;
}

DerivedRegister<Aspherity> Aspherity::reg("Aspherity");