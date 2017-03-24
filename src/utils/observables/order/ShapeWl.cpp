#include "ShapeWl.h"

ShapeWl::ShapeWl(const char* name, const char* format) : Observer(name, format) {}

ShapeWl::ShapeWl(const ShapeWl& orig) : Observer(orig) {}

ShapeWl::~ShapeWl() {}

void ShapeWl::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
    d_param = strtod(args_[ num + 1 ].c_str(), NULL);
}

double ShapeWl::observe(const Box& box, std::vector<Cell>& cells)
{
    double wlsum = 0.0;
    double N = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        wlsum += ShapeWl::calcWl(cells[i]);
        N += 1.0;
    }

    wlsum /= N;
    return wlsum;
}

double ShapeWl::calcWl(Cell& cell)
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

    double wlval = 0.0;

    if (qss > 1e-3)
    {
        wlval = wl (i_param, qlRe, qlIm) / (qss * qss * qss);
    }

    free (qlRe);
    free (qlIm);
    free (x);
    free (y);
    free (z);
    return wlval;
}

DerivedRegister<ShapeWl> ShapeWl::reg("ShapeWl");