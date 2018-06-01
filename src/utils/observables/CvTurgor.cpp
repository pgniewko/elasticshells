#include "CvTurgor.h"

CvTurgor::CvTurgor(const char* name, const char* format) : Observer(name, format) {}

CvTurgor::CvTurgor(const CvTurgor& orig) : Observer(orig) {}

CvTurgor::~CvTurgor() {}

void CvTurgor::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double CvTurgor::observe(const Box& box, std::vector<Shell>& cells, const DomainList& dl)
{
    uint N = cells.size();
    double av_turgor = 0.0;
    double av_turgor2 = 0.0;

    for (uint i = 0; i < N; i++)
    {
        av_turgor  += cells[i].getTurgor();
        av_turgor2 += cells[i].getTurgor() * cells[i].getTurgor();
    }

    av_turgor  /= N;
    av_turgor2 /= N;

    double stdev = sqrt( av_turgor2 - av_turgor * av_turgor );
    double cv = stdev / av_turgor;
    return cv;
}

DerivedRegister<CvTurgor> CvTurgor::reg("CvTurgor");