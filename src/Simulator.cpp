#include "Simulator.h"

using namespace std;

inline Vector3D Maxwell_Boltzmann(double m, const Vector3D& mean_v, double kT)
{
    return sqrt(kT / m) * Vector3D(normal(), normal(), normal()) + mean_v;
}

void compute_pairs(Pairlist& pl, Distance& domain, vector<Vector3D>& r, double RCUT, vector<Pair>& pairs)
{
    Vector3D pos[r.size()];

    for (int i = 0; i < r.size(); i++)
    {
        pos[i] = r[i];
    }


    pl.compute(domain, pos, r.size());
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
                    Vector3D r_kl = domain.delta(r[i], r[j]);
                    double R_kl = r_kl.length();

                    if (R_kl < RCUT && i != j)
                    {
                        pairs.push_back(Pair(i, j));
                    }
                }
            }
        }
}

Simulator::Simulator (const arguments& par, const Distance& _domain)
    : params(par), domain(_domain) , pl(Pairlist(par.pair_dist, 1))
{
    np = read_pos(par.input_file);
    fp.resize(np, Fields());
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
        fp[i].P = params.mass * Maxwell_Boltzmann(params.mass, Vector3D(), T);
        //cout << fp[i].P << endl;
    }

    Vector3D P;
    state(T, P);

    for (int i = 0; i < np; i++)
    {
        fp[i].P -= P;    //remove average momentum
    }

    compute_pairs(pl, domain, r, par.rcut, pairs);
}

Simulator::~Simulator()
{
}

int Simulator::read_pos(const char* filename)
{
    ifstream is(filename);
    char buf[100];

    while ( is.getline(buf, 100) )
    {
        float x, y, z;
        sscanf(buf, "%f %f %f", &x, &y, &z);
        r.push_back(Vector3D(x, y, z));

    }

    is.close();
    return r.size();
}

void Simulator::write_pos(const char* filename, bool wrap)
{
    ofstream os(filename);

    //cout << filename << " " << np<< endl;
    for (int i = 0; i < np; i++)
    {
        Vector3D pos;

        if (wrap)
        {
            pos = domain.image(r[i]);
        }
        else
        {
            pos = r[i];
        }

        os << pos.x << " " << pos.y << " " << pos.z << "\n";
    }

    os.close();
}

void Simulator::state(double& T, Vector3D& P) const
{
    double tot = 0.0;
    Vector3D tmp;

    for (int i = 0; i < np; i++)
    {
        tot += 0.5 * fp[i].P * fp[i].P / params.mass;
        tmp += fp[i].P;
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

double Simulator::compute_forces(double dt, vector<Fields>& F, int NC)
{
    double Ec = 0.0;
    Vector3D Fr_kl, Fd_kl, Fc_kl;

    for (int i = 0; i < np; i++)
    {
        F[i] = Fields();
    }

    double sqrtdt = sqrt(dt);

    for (int i = 0; i < pairs.size(); i++)
    {
        int k = pairs[i].k;
        int l = pairs[i].l;
        Vector3D r_kl = domain.delta(r[k], r[l]);
        double R_kl = r_kl.length();

        if (R_kl == 0.0)
        {
            cout << "R_KL=0" << endl;
        }

        if (R_kl < params.rcut)
        {
            //k
            Fields fk = fp[k];
            Vector3D U_k = fk.P / params.mass;
            //l
            Fields fl = fp[l];
            Vector3D U_l = fl.P / params.mass;
            //geometry
            Vector3D e_kl(r_kl / R_kl);
            Vector3D U_kl(U_k - U_l);
            double w_R = omega_R(R_kl, params.rcut);
            double w_D = w_R * w_R;

            Fc_kl = params.a * cforce(R_kl, params.rcut) * e_kl;
            Ec += 0.5 * params.a * cforce(R_kl, params.rcut) * cforce(R_kl, params.rcut);
            Fr_kl = double(NC) * params.sigma * w_R * normal() * e_kl / sqrtdt;
            Fd_kl = - double(NC) * params.gamma * w_D * (U_kl * e_kl) * e_kl;
            F[k].P +=  Fc_kl + Fd_kl + Fr_kl;
            F[l].P += -Fc_kl - Fd_kl - Fr_kl;
        }
    }

    return Ec;
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

    //else abort();

    for (int ii = start; ii <= end; ii++)
    {
        int i = abs(ii);
        int k = pairs[i].k;
        int l = pairs[i].l;
        Vector3D r_kl = domain.delta(r[k], r[l]);
        double R_kl = r_kl.length();

        if (R_kl == 0.0)
        {
            cout << "R_KL=0" << endl;
        }

        if (R_kl < params.rcut)
        {
            //k
            Fields fk = fp[k];
            Vector3D P_k = fk.P;
            //l
            Fields fl = fp[l];
            Vector3D P_l = fl.P;
            //geometry
            Vector3D e_kl(r_kl / R_kl);
            Vector3D P_kl(P_k - P_l);
            double w_R = omega_R(R_kl, params.rcut);
            double w_D = w_R * w_R;

            if (params.gamma > 0 && params.sigma > 0)
            {
                double gamma1 = 2.0 * params.gamma * w_D;
                Fr_kl = params.sigma * w_R * e_kl * sqrt((1.0 - exp(-2.0 * gamma1 * dt)) / (2.0 * gamma1) ) * normal();
                Fd_kl = 0.5 * ( P_kl * e_kl - params.a * cforce(R_kl, params.rcut) / (params.gamma * w_D) ) * e_kl * ( exp(-2.0 * params.gamma * w_D * dt) - 1.0);
            }
            else
            {
                Fr_kl *= 0.0;
                Fd_kl = params.a * cforce(R_kl, params.rcut) * e_kl * dt;
            }

            Ec += 0.5 * params.a * cforce(R_kl, params.rcut) * cforce(R_kl, params.rcut);
            fp[k].P += Fd_kl + Fr_kl;
            fp[l].P += -Fd_kl - Fr_kl;
        }
    }

    return Ec;
}

void Simulator::integrate_euler()
{
    vector<Fields> F(np, Fields());
    double dt = params.dt;
    compute_pairs(pl, domain, r, params.rcut, pairs);
    compute_forces(dt, F);

    for (int k = 0; k < np; k++)
    {
        r[k] += fp[k].P / params.mass * dt;
    }

    for (int k = 0; k < np; k++)
    {
        fp[k].P += F[k].P * dt;
    }
}

void Simulator::integrate_DPD_VV()
{
    static vector<Fields> F(np, Fields());
    vector<Vector3D> tmp_P(np, Vector3D());
    vector<Vector3D> tmp_FP(np, Vector3D());

    double dt = params.dt;
    double m = params.mass;

    for (int k = 0; k < np; k++)
    {
        r[k] += fp[k].P / m * dt + 0.5 * F[k].P / m * dt * dt;
    }

    for (int k = 0; k < np; k++)
    {
        tmp_P[k] = fp[k].P;
        tmp_FP[k] = F[k].P;
    }

    for (int k = 0; k < np; k++)
    {
        fp[k].P += 0.5 * F[k].P * dt;
    }

    compute_pairs(pl, domain, r, params.rcut, pairs);
    double Ec = compute_forces(dt, F);

    for (int k = 0; k < np; k++)
    {
        fp[k].P = tmp_P[k] + 0.5 * (F[k].P + tmp_FP[k]) * dt;
    }
}

void Simulator::integrate_trotter()
{
    double dt2 = 0.5 * params.dt;
    double dt =  params.dt;
    double m = params.mass;

    double Ec = compute_forces_trotter(dt2, 0); //2

    for (int k = 0; k < np; k++)
    {
        r[k] += fp[k].P / m * dt;    //1
    }

    compute_pairs(pl, domain, r, params.rcut, pairs);

    Ec = compute_forces_trotter(dt2, 1); //2
}