#include "AverageTurgor.h"

AverageTurgor::AverageTurgor(const char* name, const char* format) : Observer(name, format) {}

AverageTurgor::AverageTurgor(const AverageTurgor& orig) : Observer(orig) {}

AverageTurgor::~AverageTurgor() {}

double AverageTurgor::observe(const Box& box, const std::vector<Shell>& shells)
{
    uint N = shells.size();
    double av_turgor = 0.0;

    for (uint i = 0; i < N; i++)
    {
        av_turgor += shells[i].get_turgor();
    }

    av_turgor /= N;
    return av_turgor;
}

DerivedRegister<AverageTurgor> AverageTurgor::reg("AverageTurgor");