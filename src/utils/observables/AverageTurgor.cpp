#include "AverageTurgor.h"

AverageTurgor::AverageTurgor(const char* name, const char* format) : Observer(name, format) {}

AverageTurgor::AverageTurgor(const AverageTurgor& orig) : Observer(orig) {}

AverageTurgor::~AverageTurgor() {}

void AverageTurgor::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double AverageTurgor::observe(const Box& box, std::vector<Shell>& cells, const DomainList& dl)
{
    int N = cells.size();
    double av_turgor = 0.0;

    for (int i = 0; i < N; i++)
    {
        av_turgor += cells[i].getTurgor();
    }

    av_turgor /= N;
    return av_turgor;
}

DerivedRegister<AverageTurgor> AverageTurgor::reg("AverageTurgor");