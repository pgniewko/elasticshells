#include "ShapeWl.h"

ShapeWl::ShapeWl(const char* name, const char* format) : Observer(name, format) {}

ShapeWl::ShapeWl(const ShapeWl& orig) : Observer(orig) {}

ShapeWl::~ShapeWl() {}

void ShapeWl::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
}

double ShapeWl::observe(const Box& box, std::vector<Cell>& cells)
{
    double wlsum = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        wlsum += ShapeWl::calcWl(cells[i]);
    }

    wlsum /= (double) cells.size();
    return wlsum;
}

double ShapeWl::calcWl(Cell& cell)
{
    int nk = cell.getNumberVertices();
    
    double qss = 0.0;
    double* x, *y, *z;
    double* qlRe, *qlIm;
    x = (double*) malloc (nk * sizeof (double));
    y = (double*) malloc (nk * sizeof (double));
    z = (double*) malloc (nk * sizeof (double));

    for (int i = 0; i < nk ; i++)
    {
        x[i] = cell.vertices[i].r_c.x;
        y[i] = cell.vertices[i].r_c.y;
        z[i] = cell.vertices[i].r_c.z;
    }
    
    Vector3D cm = cell.getCm();
    double xc = cm.x;
    double yc = cm.y;
    double zc = cm.z;

    qlRe = (double*) malloc ((i_param + 1) * sizeof (double));
    qlIm = (double*) malloc ((i_param + 1) * sizeof (double));
    
    qlm (i_param, nk, xc, yc, zc, x, y, z, qlRe, qlIm);
    qss = qsum (i_param, qlRe, qlIm);

    double wlval = 0.0;

    if (qss > 1e-3)
    {
        wlval = Wl (i_param, qlRe, qlIm) / (qss * qss * qss);
    }

    free (qlRe);
    free (qlIm);
    free (x);
    free (y);
    free (z);
    return wlval;
}

DerivedRegister<ShapeWl> ShapeWl::reg("ShapeWl");