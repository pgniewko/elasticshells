#include "AverageVolume.h"

AverageVolume::AverageVolume(const char* name, const char* format) : Observer(name, format){}

AverageVolume::AverageVolume(const AverageVolume& orig) : Observer(orig){}

AverageVolume::~AverageVolume() {}

void AverageVolume::set_params(int num, ...)
{
    va_list arguments;
    va_start (arguments, num);
    d_param = va_arg(arguments, double);
    va_end( arguments );
};

void AverageVolume::set_params(int num, std::vector<std::string> args_)
{
    d_param = strtod(args_[ num + 0 ].c_str(), NULL);
};

double AverageVolume::observe(Box& box, std::vector<Cell>& cells)
{
    int N = cells.size();
    double av_turgor = 0.0;

    for (int i = 0; i < N; i++)
    {
        av_turgor += cells[i].calcVolume(d_param);
    }

    av_turgor /= N;
    return av_turgor;
}

DerivedRegister<AverageVolume> AverageVolume::reg("AverageVolume");

