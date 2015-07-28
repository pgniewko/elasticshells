#include "QL.h"

QL::QL(const char* name, const char* format) : Observer(name, format)
{
}

QL::QL(const QL& orig) : Observer(orig.observer_name, orig.output_format), l(orig.l), rc(orig.rc)
{
}

QL::~QL() {}

void QL::set_params(int num, ...)
{
    va_list arguments;
    va_start (arguments, num);
    l = va_arg(arguments, int);
    rc = va_arg(arguments, double);
    va_end( arguments );
}

void QL::set_params(int num, std::vector<std::string> args_)
{
    l = atoi(args_[ num+0 ].c_str());
    rc = strtod(args_[ num+1 ].c_str(), NULL);
    std::cout << "l="<< l << std::endl;
    std::cout << "rc="<< rc << std::endl;
}

double QL::observe(Box& boxs, std::vector<Cell>& cells)
{
    double qlsum = 0.0;
    double N = 0.0;

    for (int i = 0; i < cells.size(); i++)
    {
        qlsum += QL::calcQl(cells[i]);
        N += 1.0;
    }

    
    qlsum /= N;
    std::cout << qlsum<< std::endl;
    return qlsum;
}

double QL::calcQl(Cell& cell)
{
    int count;
    int n = cell.getNumberVertices();
    double qss = 0.0;
    double* x, *y, *z;
    double* qlRe, *qlIm;
    double* qlval;
    x = (double*) malloc (n * sizeof (double));
    y = (double*) malloc (n * sizeof (double));
    z = (double*) malloc (n * sizeof (double));

    for (int i = 0; i < n ; i++)
    {
        x[i] = cell.vertices[i].xyz.x;
        y[i] = cell.vertices[i].xyz.y;
        z[i] = cell.vertices[i].xyz.z;
    }

    qlval = (double*) malloc ( sizeof (double) );
    qlRe = (double*) malloc ((l + 1) * sizeof (double));
    qlIm = (double*) malloc ((l + 1) * sizeof (double));
    count = qlm (l, n, rc, x, y, z, qlRe, qlIm);
    qss = qsum (l, qlRe, qlIm);

    if (qss > 1e-3)
    {
        qlval[0] = Ql (l, count, qlRe, qlIm);
    }
    else
    {
        qlval[0] = 0.0;
    }

    free (qlRe);
    free (qlIm);
    free (x);
    free (y);
    free (z);
    return qlval[0];
}

DerivedRegister<QL> QL::reg("QL");