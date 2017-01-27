#include "QL.h"

QL::QL(const char* name, const char* format) : Observer(name, format) {}

QL::QL(const QL& orig) : Observer(orig) {}

QL::~QL() {}

void QL::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
    d_param = strtod(args_[ num + 1 ].c_str(), NULL);
}

double QL::observe(const Box& box, std::vector<Cell>& cells)
{
    double qlsum = 0.0;
    double N = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        qlsum += QL::calcQl(cells[i]);
        N += 1.0;
    }

    qlsum /= N;
    return qlsum;
}

double QL::calcQl(Cell& cell)
{
    int count;
    int n = cell.getNumberVertices();
    double qss = 0.0;
    double* x, *y, *z;
    double* qlRe, *qlIm;
    x = (double*) malloc (n * sizeof (double));
    y = (double*) malloc (n * sizeof (double));
    z = (double*) malloc (n * sizeof (double));

    for (int i = 0; i < n ; i++)
    {
        x[i] = cell.vertices[i].r_c.x;
        y[i] = cell.vertices[i].r_c.y;
        z[i] = cell.vertices[i].r_c.z;
    }

    qlRe = (double*) malloc ((i_param + 1) * sizeof (double));
    qlIm = (double*) malloc ((i_param + 1) * sizeof (double));
    double r_cutoff = cell.getInitR() * d_param;
    count = qlm (i_param, n, r_cutoff, x, y, z, qlRe, qlIm);
    qss = qsum (i_param, qlRe, qlIm);

    double qlval = 0.0;

    if (qss > 1e-3)
    {
        qlval = Ql (i_param, count, qlRe, qlIm);
    }

    free (qlRe);
    free (qlIm);
    free (x);
    free (y);
    free (z);
    return qlval;
}

DerivedRegister<QL> QL::reg("QL");