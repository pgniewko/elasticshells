#include "Integrator.h"

utils::Logger Integrator::integrator_logs("integrator");


int Integrator::FIRE_Nmin(5);
int Integrator::FIRE_N(0);
double Integrator::FIRE_DT(0.0);
double Integrator::FIRE_ALPHA(0.1);
double Integrator::FIRE_DTMAX(0.0);

Integrator::Integrator() {}

Integrator::Integrator(Simulator* s)
{
    
    integrator = &Integrator::fireIntegrator;

    FIRE_DT = s->params.dt;
    FIRE_DTMAX = 25.0 * s->params.dt;
    FIRE_ALPHA = 0.1;
    FIRE_Nmin = 5;
    FIRE_N = 0;
    integrator_logs << utils::LogLevel::FINE  << "FIRE INTEGRATOR IS USED\n";
        
    
    //setIntegrator(s);
}

Integrator::Integrator(const Integrator& orig)
{
    //TODO: implement destructor
    // delete a pointer
}

Integrator::~Integrator() {}

void Integrator::integrate(Simulator* s)
{
    (*this.*integrator)(s);
}

//void Integrator::setIntegrator(void (Integrator::*functoall)(Simulator*))
//{
//    integrator = functoall;
//}

//void Integrator::setIntegrator(Simulator* s, char* token)
//{
//    if (STRCMP (token, "fe"))
//    {
//        this->setIntegrator(&Integrator::eulerIntegrator);
//        integrator_logs << utils::LogLevel::FINE  << "FORWARD-EULER INTEGRATOR IS USED\n";
//    }
//    else if (STRCMP (token, "hm"))
//    {
//        this->setIntegrator(&Integrator::heunIntegrator);
//        integrator_logs << utils::LogLevel::FINE  << "HEUN INTEGRATOR IS USED\n";
//    }
//    else if (STRCMP (token, "rk"))
//    {
//        this->setIntegrator(&Integrator::rungeKuttaIntegrator);
//        integrator_logs << utils::LogLevel::FINE  << "RUNGE-KUTTA INTEGRATOR IS USED\n";
//    }
//    else if (STRCMP (token, "cp"))
//    {
//        this->setIntegrator(&Integrator::gearCpIntegrator);
//        integrator_logs << utils::LogLevel::FINE  << "GEAR-CP INTEGRATOR IS USED\n";
//    }
//    else if (STRCMP (token, "fire"))
//    {
//        this->setIntegrator(&Integrator::fireIntegrator);
//
//        FIRE_DT = s->params.dt;
//        FIRE_DTMAX = 25.0 * s->params.dt;
//        FIRE_ALPHA = 0.1;
//        FIRE_Nmin = 5;
//        FIRE_N = 0;
//        integrator_logs << utils::LogLevel::FINE  << "FIRE INTEGRATOR IS USED\n";
//    }
//    else
//    {
//        this->setIntegrator(&Integrator::eulerIntegrator);
//        integrator_logs << utils::LogLevel::FINE  << "DEFAULT FORWARD EULER IS USED\n";
//    }
//
//    if (integrator == NULL)
//    {
//        integrator_logs << utils::LogLevel::CRITICAL << "integrator == NULL";
//        exit(EXIT_FAILURE);
//    }
//}

void Integrator::resetParams(Simulator* s)
{
    FIRE_DT = s->params.dt;
    FIRE_ALPHA = 0.1;
    FIRE_N = 0;

    for (int i = 0; i < s->number_of_shells; i++)
    {
        for (int j = 0; j < s->shells[i].getNumberVertices(); j++)
        {
            s->shells[i].vertices[j].v_c *= 0.0; // freeze the system
            //s->shells[i].vertices[j].a_c *= 0.0; // freeze the system
        }
    }

}

/*
 * INTEGRATORS
 *
 * Viscosity of each vertex is assumed to be 1.0 !
 *
 */

void Integrator::fireIntegrator(Simulator* s)
{
    double f_inc = 1.1;
    double f_dec = 0.5;
    double a_start = 0.1;
    double f_a = 0.99;

    // MD step
    _vv(s);

    // CALC P PARAMETER
    double P = 0.0;

    for (int i = 0; i < s->number_of_shells; i++)
    {
        for (int j = 0; j < s->shells[i].getNumberVertices(); j++)
        {
            P += dot( s->shells[i].vertices[j].f_c, s->shells[i].vertices[j].v_c);
        }
    }

    //=========================

    for (int i = 0; i < s->number_of_shells; i++)
    {
        for (int j = 0; j < s->shells[i].getNumberVertices(); j++)
        {
            double v_length = s->shells[i].vertices[j].v_c.length();
            Vector3D F = s->shells[i].vertices[j].f_c;
            F.normalize();

            s->shells[i].vertices[j].v_c *= (1 - FIRE_ALPHA);
            s->shells[i].vertices[j].v_c += FIRE_ALPHA * F * v_length;
        }
    }


    if (P > 0 && FIRE_N > FIRE_Nmin)
    {
        FIRE_DT = std::min(FIRE_DTMAX, FIRE_DT * f_inc);
        FIRE_ALPHA *= f_a;
    }

    FIRE_N++;

    if (P <= 0.0)
    {
        FIRE_DT *= f_dec;
        FIRE_ALPHA = a_start;

        for (int i = 0; i < s->number_of_shells; i++)
        {
            for (int j = 0; j < s->shells[i].getNumberVertices(); j++)
            {
                s->shells[i].vertices[j].v_c *= 0.0; // freeze the system
                //s->shells[i].vertices[j].a_c *= 0.0; // freeze the system
            }
        }

        FIRE_N = 0;
    }
}

void Integrator::_vv(Simulator* s)
{
    double dt = FIRE_DT;

    // UPDATE POSITIONS
    for (int i = 0; i < s->number_of_shells; i++)
    {
        for (int j = 0; j < s->shells[i].getNumberVertices(); j++)
        {
            s->shells[i].vertices[j].r_c += dt * s->shells[i].vertices[j].v_c + 0.5 * dt * dt * s->shells[i].vertices[j].f_p; // use previously calculated forces
        }
    }

    // UPDATE VELOCITIES
    for (int i = 0; i < s->number_of_shells; i++)
    {
        for (int j = 0; j < s->shells[i].getNumberVertices(); j++)
        {
            s->shells[i].vertices[j].v_c += 0.5 * dt * s->shells[i].vertices[j].f_c;
        }
    }

    s->calcForces();

    // UPDATE VELOCITIES
    for (int i = 0; i < s->number_of_shells; i++)
    {
        for (int j = 0; j < s->shells[i].getNumberVertices(); j++)
        {
            s->shells[i].vertices[j].v_c += 0.5 * dt * s->shells[i].vertices[j].f_c;
            s->shells[i].vertices[j].f_p = s->shells[i].vertices[j].f_c; // copy forces for the next time step integration
        }
    }
}