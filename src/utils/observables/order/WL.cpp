#include "WL.h"

WL::WL(const char* name, const char* format) : Observer(name, format)
{}

WL::WL(const WL& orig) : Observer(orig.observer_name, orig.output_format), l(orig.l), rc(orig.rc)
{}

WL::~WL() {}

void WL::set_params(int num, ...)
{
    va_list arguments;
    va_start (arguments, num);
    l = va_arg(arguments, int);
    rc = va_arg(arguments, double);
    va_end( arguments );
}

void WL::set_params(int num, std::vector<std::string> args_)
{
    l = atoi(args_[ num + 0 ].c_str());
    rc = strtod(args_[ num + 1 ].c_str(), NULL);
}

double WL::observe(Box& boxs, std::vector<Cell>& cells)
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
    double* wlval;
    x = (double*) malloc (n * sizeof (double));
    y = (double*) malloc (n * sizeof (double));
    z = (double*) malloc (n * sizeof (double));

    for (int i = 0; i < n ; i++)
    {
        x[i] = cell.vertices[i].xyz.x;
        y[i] = cell.vertices[i].xyz.y;
        z[i] = cell.vertices[i].xyz.z;
    }

    wlval = (double*) malloc ( sizeof (double) );
    qlRe = (double*) malloc ((l + 1) * sizeof (double));
    qlIm = (double*) malloc ((l + 1) * sizeof (double));
    count = qlm (l, n, rc, x, y, z, qlRe, qlIm);
    qss = qsum (l, qlRe, qlIm);

    if (qss > 1e-3)
    {
        wlval[0] = wl (l, qlRe, qlIm) / (qss * qss * qss);
    }
    else
    {
        wlval[0] = 0.0;
    }

    free (qlRe);
    free (qlIm);
    free (x);
    free (y);
    free (z);
    return wlval[0];
}

DerivedRegister<WL> WL::reg("WL");