#include "Integrator.h"

utils::Logger Integrator::integrator_logs("integrator");


int Integrator::FIRE_Nmin(5);
int Integrator::FIRE_N(0);
double Integrator::FIRE_DT(0.0);
double Integrator::FIRE_ALPHA(0.1);
double Integrator::FIRE_DTMAX(0.0);

Integrator::Integrator() {}

Integrator::Integrator(Simulator* s, char* token_) 
{
    setIntegrator(s, token_);
}

Integrator::Integrator(const Integrator& orig)
{
//TODO: zaimplementuj destructor
}

Integrator::~Integrator()  {}

void Integrator::integrate(Simulator* s)
{
    (*this.*integrator)(s);
}

void Integrator::setIntegrator(void (Integrator::*functoall)(Simulator*))
{
    integrator = functoall;
}

void Integrator::resetParams(Simulator* s)
{
    FIRE_DT = s->params.dt;
    FIRE_ALPHA = 0.1;
    FIRE_N = 0;
}

void Integrator::setIntegrator(Simulator* s, char* token)
{
    if (STRCMP (token, "fe"))
    {
        this->setIntegrator(&Integrator::eulerIntegrator);
        integrator_logs << utils::LogLevel::FINE  << "FORWARD-EULER INTEGRATOR IS USED\n";
    }
    else if (STRCMP (token, "hm"))
    {
        this->setIntegrator(&Integrator::heunIntegrator);
        integrator_logs << utils::LogLevel::FINE  << "HEUN INTEGRATOR IS USED\n";
    }
    else if (STRCMP (token, "rk"))
    {
        this->setIntegrator(&Integrator::rungeKuttaIntegrator);
        integrator_logs << utils::LogLevel::FINE  << "RUNGE-KUTTA INTEGRATOR IS USED\n";
    }
    else if (STRCMP (token, "cp"))
    {
        this->setIntegrator(&Integrator::gearCpIntegrator);
        integrator_logs << utils::LogLevel::FINE  << "GEAR-CP INTEGRATOR IS USED\n";
    }
    else if (STRCMP (token, "fire"))
    {
        this->setIntegrator(&Integrator::fireIntegrator);
        
        FIRE_DT = s->params.dt;
        FIRE_DTMAX = 15.0 * s->params.dt;
        FIRE_ALPHA = 0.1;
        FIRE_Nmin = 5;
        FIRE_N = 0;
        integrator_logs << utils::LogLevel::FINE  << "FIRE INTEGRATOR IS USED\n";
    }
    else
    {
        this->setIntegrator(&Integrator::eulerIntegrator);
        integrator_logs << utils::LogLevel::FINE  << "DEFAULT FORWARD EULER IS USED\n";
    }

    if (integrator == NULL)
    {
        integrator_logs << utils::LogLevel::CRITICAL << "integrator == NULL";
        exit(EXIT_FAILURE);
    }
}


    void Integrator::eulerIntegrator(Simulator* s)
    {
        s->calcForces();
        double dt = s->params.dt;

        for (int i = 0; i < s->number_of_cells; i++)
        {
            for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
            {
                s->cells[i].vertices[j].r_c += dt * s->cells[i].vertices[j].f_c;
            }
        }
    }
        
    void Integrator::heunIntegrator(Simulator* s)
    {
    s->calcForces();
    double dt = s->params.dt;

    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            s->cells[i].vertices[j].r_p = s->cells[i].vertices[j].r_c;
            s->cells[i].vertices[j].f_p = s->cells[i].vertices[j].f_c;
        }
    }

    //move the whole time-step and calculate  forces
    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            s->cells[i].vertices[j].r_c += dt * s->cells[i].vertices[j].f_c;
        }
    }

    s->calcForces();

    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            s->cells[i].vertices[j].r_c = s->cells[i].vertices[j].r_p + 0.5 * dt * ( s->cells[i].vertices[j].f_p + s->cells[i].vertices[j].f_c);
        }
    }  
    }
        
    void Integrator::rungeKuttaIntegrator(Simulator* s)
    {
    double dt = s->params.dt;

    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            s->cells[i].vertices[j].r_p = s->cells[i].vertices[j].r_c;
        }
    }

    s->calcForces();

    //move half time-step and calculate  forces
    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            s->cells[i].vertices[j].r_c += 0.5 * dt * s->cells[i].vertices[j].f_c;
        }
    }

    s->calcForces();

    // Move the whole time-step upon the forces acting in the half-time-step
    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            s->cells[i].vertices[j].r_c = s->cells[i].vertices[j].r_p + dt * s->cells[i].vertices[j].f_c;
        }
    } 
    }
        
    void Integrator::gearCpIntegrator(Simulator* s)
    {
    double dt = s->params.dt;
    double C1, C2;

    C1 = dt;
    C2 = dt * dt / 2.0;

    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            s->cells[i].vertices[j].r_p = s->cells[i].vertices[j].r_c + C1 * s->cells[i].vertices[j].v_c + C2 * s->cells[i].vertices[j].a_c;
            s->cells[i].vertices[j].v_p = s->cells[i].vertices[j].v_c + C1 * s->cells[i].vertices[j].a_c;
            s->cells[i].vertices[j].a_p = s->cells[i].vertices[j].a_c;
        }
    }
    

    s->calcForces();

    double gear0 = 5.0 / 12.0;
    double gear2 = 1.0 / 2.0;

    double CR = gear0 * C1;
    double CA = gear2 * C1 / C2;

    Vector3D corr_v;

    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            corr_v = s->cells[i].vertices[j].f_c - s->cells[i].vertices[j].v_p; // viscosity = 1.0

            s->cells[i].vertices[j].r_c = s->cells[i].vertices[j].r_p + CR * corr_v;
            s->cells[i].vertices[j].v_c = s->cells[i].vertices[j].f_c;
            s->cells[i].vertices[j].a_c = s->cells[i].vertices[j].a_p + CA * corr_v;
        }
    }    
    }
        
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
    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            P += dot( s->cells[i].vertices[j].f_c, s->cells[i].vertices[j].v_c);
        }
    }
    
    //=========================
    
    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            double v_length = s->cells[i].vertices[j].v_c.length();
            Vector3D F = s->cells[i].vertices[j].f_c;
            F.normalize();
            
            s->cells[i].vertices[j].v_c *= (1 - FIRE_ALPHA);
            s->cells[i].vertices[j].v_c += FIRE_ALPHA * F * v_length;
        }
    }
    
    
    if (P > 0 && FIRE_N > FIRE_Nmin)
    {
        FIRE_DT = std::min(FIRE_DTMAX, FIRE_DT*f_inc);
        FIRE_ALPHA *= f_a;
    }
    
    FIRE_N++;
    
    if (P <= 0.0)
    {
        FIRE_DT *= f_dec;
        FIRE_ALPHA = a_start;
        
        for (int i = 0; i < s->number_of_cells; i++)
        {
            for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
            {
                s->cells[i].vertices[j].v_c *= 0.0; // freeze the system
                s->cells[i].vertices[j].a_c *= 0.0; // freeze the system
            }
        }
        FIRE_N = 0;
    }   
    }
    
void Integrator::_vv(Simulator* s)
{
    double dt = FIRE_DT;

    // UPDATE POSITIONS
    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            s->cells[i].vertices[j].r_c += dt * s->cells[i].vertices[j].v_c + 0.5 * dt * dt * s->cells[i].vertices[j].f_p; // use previously calculated forces
        }
    }
    
    // UPDATE VELOCITIES
    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            s->cells[i].vertices[j].v_c += 0.5 * dt * s->cells[i].vertices[j].f_c;
        }
    }
    
    s->calcForces();
    // UPDATE VELOCITIES
    for (int i = 0; i < s->number_of_cells; i++)
    {
        for (int j = 0; j < s->cells[i].getNumberVertices(); j++)
        {
            s->cells[i].vertices[j].v_c += 0.5 * dt * s->cells[i].vertices[j].f_c;
            s->cells[i].vertices[j].f_p = s->cells[i].vertices[j].f_c; // copy forces for the next time step integration
        }
    }  
}

