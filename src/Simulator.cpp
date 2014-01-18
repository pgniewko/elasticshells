#include "Simulator.h"

using namespace std;

#define STRCMP(a,b) (!strcmp(a,b))

inline Vector3D maxwell_boltzmann(double m, const Vector3D& mean_v, double kT)
{
    return sqrt(kT / m) * Vector3D(normal(), normal(), normal()) + mean_v;
}

int Simulator::initialize_pos(const char* filename, int np)
{
    if (std::string(filename).empty())
    {
        return random_pos(np);
    }
    else
    {
        return read_pos(filename);
    }
}

int Simulator::read_pos(const char* filename)
{
    ifstream is(filename);
    char buf[100];

    while ( is.getline(buf, 100) )
    {
        float x, y, z;
        sscanf(buf, "%f %f %f", &x, &y, &z);
        cells.push_back(Cell(x, y, z));
    }

    is.close();
    return cells.size();
}

int Simulator::random_pos(int np)
{
    double a = box.a();
    double b = box.b();
    double c = box.c();

    for (int i = 0; i < np; i++)
    {
        double x = uniform(0, a);
        double y = uniform(0, b);
        double z = uniform(0, c);
        cells.push_back(Cell(x, y, z));
    }
    
    return np;
}

void Simulator::write_pos(const char* filename, bool wrap)
{
    ofstream os(filename);

    for (int i = 0; i < np; i++)
    {
        Vector3D pos;

        if (wrap) /* Place all particles in the central box*/
        {
            pos = box.image(cells[i].r);
        }
        else /* Use particles absolute positions*/
        {
            pos = cells[i].r;
        }

        os << pos.x << " " << pos.y << " " << pos.z << "\n";
    }

    os.close();
}

void Simulator::write_pos_traj(const char* filename, bool wrap)
{
    ofstream os(filename, ios::app);
    os << np << "\n\n";

    for (int i = 0; i < np; i++)
    {
        Vector3D pos;

        if (wrap)
        {
            pos = box.image(cells[i].r);
        }
        else
        {
            pos = cells[i].r;
        }

        os << "H " << pos.x << " " << pos.y << " " << pos.z << "\n";
    }

    os.close();
}

Simulator::Simulator (const arguments& par, const Box& _domain)
    : params(par), box(_domain) , pl(Pairlist(par.pair_dist, 1))
{

    np = initialize_pos(par.input_file, par.n_particles);

    set_integrator(params.integrator_a);
    this->set_forces();

    cells.resize(np, Cell());
    forces.resize(np, Vector3D());

    double T;

    if (params.gamma > 0)
    {
        T = 0.5 * params.sigma * params.sigma / params.gamma;
    }
    else
    {
        T = 0;
    }

    for (int i = 0; i < np; i++)
    {
        cells[i].mass = params.mass;
        cells[i].cf = Vector3D(0, 0, 0);
        cells[i].p = params.mass * maxwell_boltzmann(params.mass, Vector3D(), T);
    }

    Vector3D P;
    state(T, P);

    for (int i = 0; i < np; i++)
    {
        cells[i].p -= P;
    }

    Pairlist::compute_pairs(pl, box, cells, par.r_cut, pairs);
}

Simulator::~Simulator()
{
    delete cForce;
}

void Simulator::state(double& T, Vector3D& P) const
{
    double tot = 0.0;
    Vector3D tmp;

    for (int i = 0; i < np; i++)
    {
        tot += 0.5 * cells[i].p * cells[i].p / cells[i].mass;
        tmp += cells[i].p;
    }

    T = 2.0 / 3.0 * tot / double(np - 1);
    P = tmp / double(np);
}


inline double omega_R(double r, double rcut)
{
    return 1.0 -  r / rcut;
}

const double SQRT3 = sqrt(3.0);

inline double weakrnd()
{
    double tmp = uniform(0.0, 1.0);

    if (tmp >= 1.0 / 3.0)
    {
        return 0.0;
    }
    else if (tmp < 1.0 / 6.0)
    {
        return SQRT3;
    }
    else
    {
        return -SQRT3;
    }
}

double Simulator::compute_forces(double dt)
{
    
    double Ec = 0.0;
    Vector3D Fr_kl, Fd_kl, Fc_kl;

    for (int i = 0; i < np; i++)
    {
        forces[i] = Vector3D();
        cells[i].reset_cforce();
    }

    double sqrtdt = sqrt(dt);

    for (int i = 0; i < pairs.size(); i++)
    {
        int k = pairs[i].k;
        int l = pairs[i].l;
        Vector3D r_kl = box.delta(cells[k].r, cells[l].r);
        double R_kl = r_kl.length();

        if (R_kl == 0.0)
        {
            cout << "R_KL=0" << endl;
        }

        if (R_kl < params.r_cut)
        {
            //k
            Vector3D U_k = cells[k].p / cells[k].mass ;
            //l
            Vector3D U_l = cells[l].p / cells[l].mass;
            //geometry
            Vector3D e_kl(r_kl / R_kl);
            Vector3D U_kl(U_k - U_l);
            double w_R = omega_R(R_kl, params.r_cut);
            double w_D = w_R * w_R;

            Fc_kl = cForce->eval(r_kl, cells[k], cells[l]);   
            Fr_kl = params.sigma * w_R * normal() * e_kl / sqrtdt;
            Fd_kl = -params.gamma * w_D * (U_kl * e_kl) * e_kl;
            forces[k] +=  Fc_kl + Fd_kl + Fr_kl;
            forces[l] += -Fc_kl - Fd_kl - Fr_kl;
            
            cells[k].cf +=  Fc_kl;
            cells[l].cf += -Fc_kl;
        }
    }

    return Ec;
}

void Simulator::integrate_euler()
{
    double dt = params.dt;
    Pairlist::compute_pairs(pl, box, cells, params.r_cut, pairs);
    compute_forces(dt);

    for (int k = 0; k < np; k++)
    {
        cells[k].r += cells[k].p / cells[k].mass * dt;
    }

    for (int k = 0; k < np; k++)
    {
       cells[k].p += forces[k] * dt;
    }
}

void Simulator::integrate_DPD_VV() 
/* agreement with original mydpd - without "-ffast-math -Wno-deprecated" flags */
/* Otherwise round-off errors*/
{
    
    vector<Vector3D> tmp_P(np, Vector3D());
    vector<Vector3D> tmp_FP(np, Vector3D());

    double dt = params.dt;

    for (int k = 0; k < np; k++)
    {
        cells[k].r += cells[k].p / cells[k].mass * dt + 0.5 * forces[k] / cells[k].mass * dt * dt;
    }

    for (int k = 0; k < np; k++)
    {
        tmp_P[k] = cells[k].p;
        tmp_FP[k] = forces[k];
    }

    for (int k = 0; k < np; k++)
    {
        cells[k].p += 0.5 * forces[k] * dt;
    }

    Pairlist::compute_pairs(pl, box, cells, params.r_cut, pairs);
    double Ec = compute_forces(dt);

    for (int k = 0; k < np; k++)
    {
        cells[k].p = tmp_P[k] + 0.5 * (forces[k] + tmp_FP[k]) * dt;
    }
}


double Simulator::compute_forces_trotter(double dt, int order)
{
    double Ec = 0.0;
    Vector3D Fr_kl, Fd_kl, Fc_kl;
    double sqrtdt = sqrt(dt);

    for (int i = 0; i < np; i++)
    {
        cells[i].reset_cforce();
    }
    
    int start, end, step;

    if (order == 0)
    {
        start = 0;
        end = pairs.size() - 1;
    }
    else if (order == 1)
    {
        start = - pairs.size() + 1;
        end = 0;
    }
    else
    {
        start = 0;
        end = pairs.size() - 1;
    }

//    else abort(); /* "stdlib.h" */

    for (int ii = start; ii <= end; ii++)
    {
        int i = abs(ii);
        int k = pairs[i].k;
        int l = pairs[i].l;
        Vector3D r_kl = box.delta(cells[k].r, cells[l].r);
        double R_kl = r_kl.length();

        if (R_kl == 0.0)
        {
            cout << "R_KL=0" << endl;
        }

        if (R_kl < params.r_cut)
        {
            //k
            Vector3D P_k = cells[k].p;
            //l
            Vector3D P_l = cells[l].p;
            //geometry
            Vector3D e_kl(r_kl / R_kl);
            Vector3D P_kl(P_k - P_l);
            double w_R = omega_R(R_kl, params.r_cut);
            double w_D = w_R * w_R;

            Fc_kl = cForce->eval(r_kl, cells[k], cells[l]);
            if (params.gamma > 0 && params.sigma > 0)
            {
                double gamma1 = 2.0 * params.gamma * w_D;
                Fr_kl = params.sigma * w_R * e_kl * sqrt((1.0 - exp(-2.0 * gamma1 * dt)) / (2.0 * gamma1) ) * normal();
                Fd_kl = 0.5 * ( P_kl * e_kl - cForce->magn(r_kl, cells[k], cells[l]) / (params.gamma * w_D) ) * e_kl * ( exp(-2.0 * params.gamma * w_D * dt) - 1.0);
                
            }
            else
            {
                Fr_kl *= 0.0;
                Fd_kl = Fc_kl * dt;
            }

            //Ec += 0.5 * params.a * cforce(R_kl, params.r_cut) * cforce(R_kl, params.r_cut);
            cells[k].p +=  Fd_kl + Fr_kl;
            cells[l].p += -Fd_kl - Fr_kl;
            
            cells[k].cf +=  Fc_kl;
            cells[l].cf += -Fc_kl;
            
        }
    }

    return Ec;
}

void Simulator::integrate_trotter()
{
    double dt2 = 0.5 * params.dt;
    double dt =  params.dt;

    double Ec = compute_forces_trotter(dt2, 0); //2 

    for (int k = 0; k < np; k++)
    {
        cells[k].r += cells[k].p / cells[k].mass * dt;    //1
    }

    Pairlist::compute_pairs(pl, box, cells, params.r_cut, pairs);
    
    Ec = compute_forces_trotter(dt2, 1); //2
}

void Simulator::set_integrator(void (Simulator::*functoall)())
{
    integrator = functoall;
}

void Simulator::set_integrator(char* token)
{
    if (STRCMP (token, "vv"))
    {
        this->set_integrator(&Simulator::integrate_DPD_VV);
    }
    else if (STRCMP (token, "trot"))
    {
        this->set_integrator(&Simulator::integrate_trotter);
    }
    else
    {
        this->set_integrator(&Simulator::integrate_DPD_VV);
    }
}

void Simulator::set_forces()
{
    cForce = new RepulsiveForce(params.a, params.r_cut);
}

void Simulator::integrate()
{
    (*this.*integrator)();
}