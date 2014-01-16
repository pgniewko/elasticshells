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

    //error (123, 0, "ABORTED");
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

    //cout << filename << " " << np<< endl;
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

void compute_pairs(Pairlist& pl, Box& domain, vector<Cell>& cells, double RCUT, vector<Pair>& pairs)
{
    Vector3D pos[cells.size()];

    for (int i = 0; i < cells.size(); i++)
    {
        pos[i] = cells[i].r;
    }


    //pl.compute(domain, pos, r.size());
    pl.compute(domain, pos, cells.size());
    pairs.clear();

    for (int k = 0; k < pl.boxlist.size(); k++)//iterate over boxes
        for (int nbr = 0; nbr < pl.nbrlist[k].size(); nbr++)
        {
            //iterate over nbr boxes
            int k1 =  pl.nbrlist[k][nbr];

            for (int l = 0; l < pl.boxlist[ k ].size(); l++)  //iterate over atoms
            {
                int i =  pl.boxlist[k][l];
                int startl1 = 0;

                if (k1 == k)
                {
                    startl1 = l + 1;
                }

                for (int l1 = startl1; l1 < pl.boxlist[k1].size(); l1++)
                {
                    //iterate over atoms
                    int j = pl.boxlist[k1][l1];
                    Vector3D r_kl = domain.delta(cells[i].r, cells[j].r);
                    double R_kl = r_kl.length();

                    if (R_kl < RCUT && i != j)
                    {
                        pairs.push_back(Pair(i, j));
                    }
                }
            }
        }
}

Simulator::Simulator (const arguments& par, const Box& _domain)
    : params(par), box(_domain) , pl(Pairlist(par.pair_dist, 1))
{

    np = initialize_pos(par.input_file, par.n_particles);

    set_integrator(params.integrator_a);

    cells.resize(np, Cell());


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
        cells[i].f = Vector3D(0, 0, 0);
        //fp[i].P = params.mass * maxwell_boltzmann(params.mass, Vector3D(), T);
        cells[i].p = params.mass * maxwell_boltzmann(params.mass, Vector3D(), T);

        //cout << fp[i].P << endl;
    }

    Vector3D P;
    state(T, P);

    for (int i = 0; i < np; i++)
    {
        //fp[i].P -= P;    //remove average momentum
        cells[i].p -= P;
    }

    //compute_pairs(pl, box, r, par.r_cut, pairs, cells);
    compute_pairs(pl, box, cells, par.r_cut, pairs);
}

//Simulator::~Simulator()
//{
//}

void Simulator::state(double& T, Vector3D& P) const
{
    double tot = 0.0;
    Vector3D tmp;

    for (int i = 0; i < np; i++)
    {
        //tot += 0.5 * fp[i].P * fp[i].P / params.mass;
        //tmp += fp[i].P;

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

inline double cforce(double r, double rcut)
{
    return (1.0 -  r / rcut);
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

inline double rho_p(int p)
{
    double r = 0.0;

    for (int i = 1; i <= p; i++)
    {
        r += 1.0 / double(i * i);
    }

    //return 1.0 / 12.0 - 1.0 / (2.0 * M_PI * M_PI) * r;
    return 1.0 / 12.0 - 1.0 / (2.0 * PI * PI) * r;
}

double Simulator::compute_forces(double dt, vector<Vector3D>& F, int NC)
{
    double Ec = 0.0;
    Vector3D Fr_kl, Fd_kl, Fc_kl;

    for (int i = 0; i < np; i++)
    {
        F[i] = Vector3D();
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
            //Fields fk = fp[k];
            Vector3D U_k = cells[k].p / cells[k].mass ;//fk.P / params.mass;
            //l
            //Fields fl = fp[l];
            Vector3D U_l = cells[l].p / cells[l].mass; //fl.P / params.mass;
            //geometry
            Vector3D e_kl(r_kl / R_kl);
            Vector3D U_kl(U_k - U_l);
            double w_R = omega_R(R_kl, params.r_cut);
            double w_D = w_R * w_R;

            Fc_kl = params.a * cforce(R_kl, params.r_cut) * e_kl;
            Ec += 0.5 * params.a * cforce(R_kl, params.r_cut) * cforce(R_kl, params.r_cut);
            Fr_kl = double(NC) * params.sigma * w_R * normal() * e_kl / sqrtdt;
            Fd_kl = - double(NC) * params.gamma * w_D * (U_kl * e_kl) * e_kl;
            F[k] +=  Fc_kl + Fd_kl + Fr_kl;
            F[l] += -Fc_kl - Fd_kl - Fr_kl;
        }
    }

    return Ec;
}

void Simulator::integrate_euler()
{
    vector<Vector3D> F(np, Vector3D());
    double dt = params.dt;
    compute_pairs(pl, box, cells, params.r_cut, pairs);
    compute_forces(dt, F);

    for (int k = 0; k < np; k++)
    {
        cells[k].r += cells[k].p / cells[k].mass * dt;
    }

    for (int k = 0; k < np; k++)
    {
        cells[k].p += F[k] * dt;
    }
}

void Simulator::integrate_DPD_VV()
{
    static vector<Vector3D> F(np, Vector3D());
    vector<Vector3D> tmp_P(np, Vector3D());
    vector<Vector3D> tmp_FP(np, Vector3D());

    double dt = params.dt;
    //double m = params.mass;

    for (int k = 0; k < np; k++)
    {
        cells[k].r += cells[k].p / cells[k].mass * dt + 0.5 * F[k] / cells[k].mass * dt * dt;
    }

    for (int k = 0; k < np; k++)
    {
        tmp_P[k] = cells[k].p; // fp[k].P;
        tmp_FP[k] = F[k];
    }

    for (int k = 0; k < np; k++)
    {
        cells[k].p += 0.5 * F[k] * dt;
    }

    compute_pairs(pl, box, cells, params.r_cut, pairs);
    double Ec = compute_forces(dt, F);

    for (int k = 0; k < np; k++)
    {
        cells[k].p = tmp_P[k] + 0.5 * (F[k] + tmp_FP[k]) * dt;
    }
}


double Simulator::compute_forces_trotter(double dt, int order)
{
//	static double rhop =  rho_p(20);
    double Ec = 0.0;
    Vector3D Fr_kl, Fd_kl, Fc_kl;
    double sqrtdt = sqrt(dt);

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
            //Fields fk = fp[k];
            Vector3D P_k = cells[k].p; //fk.P;
            //l
            //Fields fl = fp[l];
            Vector3D P_l = cells[l].p; //fl.P;
            //geometry
            Vector3D e_kl(r_kl / R_kl);
            Vector3D P_kl(P_k - P_l);
            double w_R = omega_R(R_kl, params.r_cut);
            double w_D = w_R * w_R;

            if (params.gamma > 0 && params.sigma > 0)
            {
                double gamma1 = 2.0 * params.gamma * w_D;
                Fr_kl = params.sigma * w_R * e_kl * sqrt((1.0 - exp(-2.0 * gamma1 * dt)) / (2.0 * gamma1) ) * normal();
                Fd_kl = 0.5 * ( P_kl * e_kl - params.a * cforce(R_kl, params.r_cut) / (params.gamma * w_D) ) * e_kl * ( exp(-2.0 * params.gamma * w_D * dt) - 1.0);
            }
            else
            {
                Fr_kl *= 0.0;
                Fd_kl = params.a * cforce(R_kl, params.r_cut) * e_kl * dt;
            }

            Ec += 0.5 * params.a * cforce(R_kl, params.r_cut) * cforce(R_kl, params.r_cut);
            //fp[k].P += Fd_kl + Fr_kl;
            //fp[l].P += -Fd_kl - Fr_kl;
            cells[k].p += Fd_kl + Fr_kl;
            cells[l].p += -Fd_kl - Fr_kl;
        }
    }

    return Ec;
}

void Simulator::integrate_trotter()
{
    double dt2 = 0.5 * params.dt;
    double dt =  params.dt;
    //double m = params.mass;

    double Ec = compute_forces_trotter(dt2, 0); //2

    for (int k = 0; k < np; k++)
    {
        cells[k].r += cells[k].p / cells[k].mass * dt;    //1
    }

    compute_pairs(pl, box, cells, params.r_cut, pairs);

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
        this->set_integrator(&Simulator::integrate_trotter);
    }
}

void Simulator::integrate()
{
    (*this.*integrator)();
}