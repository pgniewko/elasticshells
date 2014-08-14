#include "Simulator.h"

Simulator::Simulator() {
}

Simulator::Simulator(const Simulator& orig) {
}

Simulator::~Simulator() {
}


void Simulator::setIntegrator(char* token)
{
    if (STRCMP (token, "vv"))
    {
        this->setIntegrator(&Simulator::integrate_VV);
    }
    else if (STRCMP (token, "eu"))
    {
        this->setIntegrator(&Simulator::integrate_euler);
    }
    else if (STRCMP (token, 'de'))
    {
        this->setIntegrator(&Simulator::integrate_dempedeuler);
    }
    else
    {
        this->setIntegrator(&Simulator::integrate_euler);
    }
}


void Simulator::integrate()
{
    (*this.*integrator)();
}
