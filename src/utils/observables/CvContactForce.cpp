#include "CvContactForce.h"

CvContactForce::CvContactForce(const char* name, const char* format) : Observer(name, format) {}

CvContactForce::CvContactForce(const CvContactForce& orig) : Observer(orig) {}

CvContactForce::~CvContactForce() {}

void CvContactForce::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
};

double CvContactForce::observe(const Box& box, std::vector<Cell>& cells, const DomainList& dl)
{
    double average_cf  = 0.0;
    double average_cf2 = 0.0;
    double counter = 1.0;

    uint N = cells.size();
    double cf = 0.0;

    for (uint i = 0; i < N; i++)
    {
        for (uint j = i + 1; j < N; j++)
        {
            if (i_param == 0)
            {
                cf = cells[i].contactForce(cells[j], box, true);
            }
            else if (i_param == 1)
            {
                cf = dl.calcContactForce(i, j, cells, box);
            }

            if (cf > 0.0)
            {
                average_cf  += cf;
                average_cf2 += (cf * cf);
                counter += 1.0;
            }
        }
    }

    if (counter == 0)
    {
        return 0.0;
    }

    average_cf  /= counter;
    average_cf2 /= counter;

    double stdev = sqrt( average_cf2 - average_cf * average_cf );

    if (stdev == 0.0)
    {
        return 0.0;
    }

    double cv = stdev / average_cf;
    return cv;
}

DerivedRegister<CvContactForce> CvContactForce::reg("CvContactForce");