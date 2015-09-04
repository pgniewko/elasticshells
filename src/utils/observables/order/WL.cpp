#include "WL.h"

WL::WL(const char* name, const char* format) : Observer(name, format) {}

WL::WL(const WL& orig) : Observer(orig) {}

WL::~WL() {}

void WL::set_params(int num, ...)
{
    va_list arguments;
    va_start (arguments, num);
    i_param = va_arg(arguments, int);
    d_param = va_arg(arguments, double);
    va_end( arguments );
}

void WL::set_params(int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
    d_param = strtod(args_[ num + 1 ].c_str(), NULL);
}

double WL::observe(Box& box, std::vector<Cell>& cells)
{
    double wlsum = 0.0;
    double N = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        wlsum += WL::calcWl(cells[i]);
        N += 1.0;
    }

    wlsum /= N;
    return wlsum;
}

double WL::calcWl(Cell& cell)
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
        x[i] = cell.vertices[i].xyz.x;
        y[i] = cell.vertices[i].xyz.y;
        z[i] = cell.vertices[i].xyz.z;
    }

    qlRe = (double*) malloc ((i_param + 1) * sizeof (double));
    qlIm = (double*) malloc ((i_param + 1) * sizeof (double));
    count = qlm (i_param, n, d_param, x, y, z, qlRe, qlIm);
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

DerivedRegister<WL> WL::reg("WL");