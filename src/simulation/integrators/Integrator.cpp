#include "Integrator.h"

utils::Logger Integrator::integrator_logs("integrator");


int Integrator::FIRE_Nmin(5);
int Integrator::FIRE_N(0);
int Integrator::ITER_NUM(0);
int Integrator::TOTAL_ITER_NUM(0);
double Integrator::FIRE_DT(0.0);
double Integrator::FIRE_ALPHA(0.1);
double Integrator::FIRE_DTMAX(0.0);

Integrator::Integrator() : n(0) {}

Integrator::Integrator(Simulator* s) : n(0)
{
    integrator = &Integrator::fire_integrator;

    FIRE_DT = s->params.dt;
    FIRE_DTMAX = 25.0 * s->params.dt;
    FIRE_ALPHA = 0.1;
    FIRE_Nmin = 5;
    FIRE_N = 0;
    integrator_logs << utils::LogLevel::FINE  << "FIRE INTEGRATOR IS USED\n";
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

void Integrator::reset_params(Simulator* s)
{
    FIRE_DT = s->params.dt;
    FIRE_ALPHA = 0.1;
    FIRE_N = 0;

    for (uint i = 0; i < n; i++)
    {
        vel[i] = 0.0;
    }

    for (uint i = 0; i < n; i++)
    {
        f_p[i] = s->forces[i];
    }

}

void Integrator::fire_integrator(Simulator* s)
{
    double f_inc = 1.1;
    double f_dec = 0.5;
    double a_start = 0.1;
    double f_a = 0.99;

    double fv = 0.0;
    double vv = 0.0;
    double ff = 0.0;

    double fx, fy, fz;
    for (uint i = 0; i < n; i++)
    {
        fv += s->forces[i] * vel[i];
        vv += vel[i] * vel[i];
        ff += s->forces[i] * s->forces[i];
    }

    for (uint i = 0; i < n; i += 3)
    {
        fx = s->forces[i + 0];
        fy = s->forces[i + 1];
        fz = s->forces[i + 2];
        vel[i + 0] *= (1 - FIRE_ALPHA);
        vel[i + 1] *= (1 - FIRE_ALPHA);
        vel[i + 2] *= (1 - FIRE_ALPHA);
        vel[i + 0] += FIRE_ALPHA * sqrt(vv) * fx / sqrt(ff);
        vel[i + 1] += FIRE_ALPHA * sqrt(vv) * fy / sqrt(ff);
        vel[i + 2] += FIRE_ALPHA * sqrt(vv) * fz / sqrt(ff);
    }

    if (fv > 0 && FIRE_N > FIRE_Nmin)
    {
        FIRE_DT = std::min(FIRE_DTMAX, FIRE_DT * f_inc);
        FIRE_ALPHA *= f_a;
    }

    FIRE_N++;

    if (fv <= 0.0)
    {
        FIRE_DT *= f_dec;
        FIRE_ALPHA = a_start;

        for (uint i = 0; i < n; i++)
        {
            vel[i] = 00.;
        }

        FIRE_N = 0;
    }
    
    // MD step
    _vv(s);
}

void Integrator::_vv(Simulator* s)
{
    double dt = FIRE_DT;

    // UPDATE POSITIONS
    for (uint i = 0; i < n; i++)
    {
        s->xyz[i] += dt * vel[i] + 0.5 * dt * dt * f_p[i]; // use previously calculated forces
    }

    // UPDATE VELOCITIES
    for (uint i = 0; i < n; i++)
    {
        vel[i] += 0.5 * dt * s->forces[i];
    }

    s->calculate_forces();

    // UPDATE VELOCITIES
    for (uint i = 0; i < n; i++)
    {
        vel[i] += 0.5 * dt * s->forces[i];
        f_p[i] = s->forces[i];
    }
}

void Integrator::set_n(uint n_)
{
    n = n_;

    for (uint i = 0; i < n; i++)
    {
        vel.push_back(0.0);
        f_p.push_back(0.0);
    }
}