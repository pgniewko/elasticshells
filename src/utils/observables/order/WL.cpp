#include "WL.h"

WL::WL() {}

WL::WL(const WL& orig) {}

WL::~WL() {}

double WL::calcWl(std::vector<Cell>& cells, int l, double rc)
{
    double wlsum = 0.0;
    double N = 0.0;

    for (int i = 0; i < cells.size(); i++)
    {
        wlsum += WL::calcWl(cells[i], l, rc);
        N += 1.0;
    }

    wlsum /= N;
    return wlsum;
}

double WL::calcWl(Cell& cell, int l, double rc)
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