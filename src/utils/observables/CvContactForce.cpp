#include "CvContactForce.h"

CvContactForce::CvContactForce() {
}

CvContactForce::CvContactForce(const CvContactForce& orig) {
}

CvContactForce::~CvContactForce() {
}


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
                counter + 1.0;
            }
        } 
    }
    
    average_cf  /= N;
    average_cf2 /= N;
    
    double stdev = sqrt( average_cf2 - average_cf * average_cf );
    double cv = stdev / average_cf;
    return cv;
    
    
}
