#include "AverageTurgor.h"

AverageTurgor::AverageTurgor() {}

AverageTurgor::AverageTurgor(const AverageTurgor& orig) {}

AverageTurgor::~AverageTurgor() {}

double AverageTurgor::populationAverageTurgor(std::vector<Cell>& cells)
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