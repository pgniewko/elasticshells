#include "Tensorial.h"

Tensorial::Tensorial(const char* name, const char* format) : Observer(name, format) {}

Tensorial::Tensorial(const Tensorial& orig) : Observer(orig) {}

Tensorial::~Tensorial() {}

void Tensorial::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
};

double Tensorial::observe(const Box& box, const std::vector<Shell>& shells)
{
    double Ob = 0.0;

    for (uint i = 0; i < shells.size(); i++)
    {
        double A[3][3];
        double V[3][3];
        double d[3];

        for (int i_ = 0; i_ < 3; i_++)
        {
            d[i_] = 0.0;

            for (int j_ = 0; j_ < 3; j_++)
            {
                A[i_][j_] = 0.0;
                V[i_][j_] = 0.0;
            }
        }

        Vector3D cell_cm = shells[i].center_of_mass;

        Vector3D rj;

        double ev1 = 0.0, ev2 = 0.0, ev3 = 0.0;

        for (int j = 0; j < shells[i].get_number_vertices(); j++)
        {
            rj = shells[i].vertices[j].r_c - cell_cm;

            A[0][0] += rj.x * rj.x;
            A[1][1] += rj.y * rj.y;
            A[2][2] += rj.z * rj.z;

            A[0][1] += rj.x * rj.y;
            A[0][2] += rj.x * rj.z;

            A[1][0] += rj.y * rj.x;
            A[1][2] += rj.y * rj.z;

            A[2][0] += rj.z * rj.x;
            A[2][1] += rj.z * rj.y;
        }

        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
            {
                A[j][k] /= shells[i].get_number_vertices();
            }

        eigen_decomposition(A, V, d);

        ev1 = std::max(d[0], d[1]);
        ev1 = std::max(ev1, d[2]);

        ev3 = std::min(d[0], d[1]);
        ev3 = std::min(ev3, d[2]);

        ev2 = d[0] + d[1] + d[2] - ev1 - ev3;

        if (i_param == 0) // RADIUS OF GYRATION
        {
            Ob += (ev1 + ev2 + ev3);
        }
        else if (i_param == 1) // ASPHERICITY
        {
            Ob += ( ev1 - 0.5 * (ev2 + ev3) );
        }
        else if (i_param == 2) // ACYLIDRICITY
        {
            Ob += ( ev2 - ev3 );
        }
        else if (i_param == 3) // ASPECT-RATIO
        {
            if (ev3 > 0)
            {
                Ob += (ev1 / ev3);
            }
            else
            {
                throw DataException("Eigenvalue < 0 (or equal to 0)");
            }
        }
        else
        {
            Ob += 0.0;
        }
    }

    Ob /= shells.size();

    return Ob;
}

DerivedRegister<Tensorial> Tensorial::reg("Tensorial");