#include "CvContactForce.h"

CvContactForce::CvContactForce(const char* name, const char* format) : Observer(name, format) {}

CvContactForce::CvContactForce(const CvContactForce& orig) : Observer(orig) {}

CvContactForce::~CvContactForce() {}

void CvContactForce::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double CvContactForce::observe(const Box& box, std::vector<Cell>& cells)
{
    double average_cf  = 0.0;
    double average_cf2 = 0.0;
    double counter = 1.0;
    
    uint N = cells.size();
    double cf = 0.0;
    
    for (uint i = 0; i < N; i++)
    {
        for (uint j = i+1; j < N; j++)
        {
            cf = cells[i].contactForce(cells[j], box);
            if (cf > 0.0)
            {
                average_cf  += cf;
                average_cf2 += (cf*cf);
                counter += 1.0;
            }
        } 
    }
    
    average_cf  /= counter;
    average_cf2 /= counter;
    
    double stdev = sqrt( average_cf2 - average_cf * average_cf );
    double cv = stdev / average_cf;
    return cv;
}

DerivedRegister<CvContactForce> CvContactForce::reg("CvContactForce");