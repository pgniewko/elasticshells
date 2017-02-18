#include "AverageContactForce.h"

AverageContactForce::AverageContactForce(const char* name, const char* format) : Observer(name, format) {}

AverageContactForce::AverageContactForce(const AverageContactForce& orig) : Observer(orig) {}

AverageContactForce::~AverageContactForce() {
}

void AverageContactForce::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double AverageContactForce::observe(const Box& box, std::vector<Cell>& cells)
{
    double average_cf = 0.0;
    double counter = 1.0;
    
    uint N = cells.size();
    double cf = 0.0;
    
    for (uint i = 0; i < N; i++)
    {
        for (uint j = i+1; j < N; j++)
        {
            cf = cells[i].;
            if (cf > 0.0)
            {
                cells[i].contactForce(cells[j], box);
                counter += 1.0;
            }
        } 
    }
    
    average_cf /= counter;
    return average_cf;
    
    
}


DerivedRegister<AverageContactForce> AverageContactForce::reg("AverageContactForce");