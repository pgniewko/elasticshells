#include "PackingWl.h"

PackingWl::PackingWl(const char* name, const char* format) : Observer(name, format) {}

PackingWl::PackingWl(const PackingWl& orig) : Observer(orig) {}

PackingWl::~PackingWl() {}

void PackingWl::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
}

double PackingWl::observe(const Box& box, std::vector<Cell>& cells, const DomainList& dl)
{
   if (box.pbc == false)
    {
        return 0.0;    
    }
    
    double qlsum = 0.0;

    for (uint i = 0; i < cells.size(); i++)
    {
        qlsum += PackingWl::calcWl(box, cells, i, dl);
    }

    qlsum /= (double) cells.size();
    return qlsum;
}

double PackingWl::calcWl(const Box& box, std::vector<Cell>& cells, unsigned int ci, const DomainList& dl)
{
    Vector3D ci_cm = cells[ci].getCm();
    std::vector<Vector3D> neighs_cm;
    for (uint i = 0; i < cells.size(); i++)
    {
        if (i != ci)
        {
            //if ( cells[ci].isInContact(cells[i], box) )
            if ( dl.isInContact(ci, i, cells, box) );
            {
                Vector3D dij;
                Vector3D cm = cells[i].getCm();
                Box::getDistance(dij, ci_cm, cm, box);
                cm = ci_cm + dij;
                neighs_cm.push_back(cm);
            }
        }
    }
    
    int nk = neighs_cm.size();
    
    double qss = 0.0;
    double *x, *y, *z;
    double *qlRe, *qlIm;
    x = (double*) malloc (nk * sizeof (double));
    y = (double*) malloc (nk * sizeof (double));
    z = (double*) malloc (nk * sizeof (double));

    for (int i = 0; i < nk ; i++)
    {
        x[i] = neighs_cm[i].x;
        y[i] = neighs_cm[i].y;
        z[i] = neighs_cm[i].z;
    }

    double xc = ci_cm.x;
    double yc = ci_cm.y;
    double zc = ci_cm.z;

    qlRe = (double*) malloc ((i_param + 1) * sizeof (double));
    qlIm = (double*) malloc ((i_param + 1) * sizeof (double));
    
    qlm(i_param, nk, xc, yc, zc, x, y, z, qlRe, qlIm);
    qss = qsum(i_param, qlRe, qlIm);

    double wlval = 0.0;

    if (qss > 1e-4)
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

DerivedRegister<PackingWl> PackingWl::reg("PackingWl");