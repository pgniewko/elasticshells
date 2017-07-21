#include "ShapeQl.h"

ShapeQl::ShapeQl(const char* name, const char* format) : Observer(name, format) {}

ShapeQl::ShapeQl(const ShapeQl& orig) : Observer(orig) {}

ShapeQl::~ShapeQl() {}

void ShapeQl::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
}

double ShapeQl::observe(const Box& box, std::vector<Cell>& cells, const DomainList& dl)
{
    double qlsum = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        qlsum += ShapeQl::calcQl(cells[i]);
    }

    qlsum /= (double) cells.size();
    return qlsum;
}

double ShapeQl::calcQl(Cell& cell)
{
    int nk = cell.getNumberVertices();

    if (nk == 1)
    {
        return 0.0;
    }

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

    qlm(i_param, nk, xc, yc, zc, x, y, z, qlRe, qlIm);
    qss = qsum(i_param, qlRe, qlIm);

    double qlval = 0.0;

    if (qss > 1e-4)
    {
        qlval = Ql(i_param, nk, qlRe, qlIm);
    }

    free (qlRe);
    free (qlIm);
    free (x);
    free (y);
    free (z);
    return qlval;
}

DerivedRegister<ShapeQl> ShapeQl::reg("ShapeQl");