#include "PackingQl.h"

PackingQl::PackingQl(const char* name, const char* format) : Observer(name, format) {}

PackingQl::PackingQl(const PackingQl& orig) : Observer(orig) {}

PackingQl::~PackingQl() {}

void PackingQl::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
}

double PackingQl::observe(const Box& box, std::vector<Cell>& cells)
{
    if (box.pbc == false)
    {
        return 0.0;    
    }
    
    double qlsum = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        qlsum += PackingQl::calcQl(cells[i]);
    }

    qlsum /= (double) cells.size();
    return qlsum;
}

double PackingQl::calcQl(const Box& box, std::vector<Cell>& cells, unsigned int ci)
{
    std::vector<Vector3D> neighs_cm;
    for (uint i = 0; i < cells.size(); i++)
    {
        if (i != ci)
        {
            if ( cells[ci].isInContact(cells[i], box) )
            {
                Vector3D cm = cells[i].getCm();
                neighs_cm.push_back(cm);
            }
        }
    }
    
    int nk = neighs_cm.size();
    
    double qss = 0.0;
    double* x, *y, *z;
    double* qlRe, *qlIm;
    x = (double*) malloc (nk * sizeof (double));
    y = (double*) malloc (nk * sizeof (double));
    z = (double*) malloc (nk * sizeof (double));

    for (int i = 0; i < nk ; i++)
    {
        x[i] = neighs_cm[i].x;
        y[i] = neighs_cm[i].y;
        z[i] = neighs_cm[i].z;
    }

    Vector3D cm = cells[ci].getCm();
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

DerivedRegister<PackingQl> PackingQl::reg("PackingQl");