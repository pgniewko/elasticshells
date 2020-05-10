#include "PackingWl.h"

PackingWl::PackingWl(const char* name, const char* format) : Observer(name, format) {}

PackingWl::PackingWl(const PackingWl& orig) : Observer(orig) {}

PackingWl::~PackingWl() {}

double PackingWl::observe(const Box& box, const std::vector<Shell>& shells)
{
    if (box.pbc == false)
    {
        return -1;
    }

    double qlsum = 0.0;

    create_shells_image(box, shells);
    copy_shells_data(box, shells);

    for (uint i = 0; i < shells.size(); i++)
    {
        qlsum += PackingWl::calc_wl(box, shells, i);
    }

    qlsum /= (double) shells.size();
    return qlsum;
}

double PackingWl::calc_wl(const Box& box, const std::vector<Shell>& shells, unsigned int ci)
{
    Vector3D ci_cm = shells[ci].get_cm();
    std::vector<Vector3D> neighs_cm;

    for (uint i = 0; i < shells.size(); i++)
    {
        if (i != ci)
        {
            if ( is_in_contact(i, ci) )
            {
                Vector3D dij;
                Vector3D cm = shells[i].get_cm();
                Box::get_distance(dij, ci_cm, cm, box);
                cm = ci_cm + dij;
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

    double xc = ci_cm.x;
    double yc = ci_cm.y;
    double zc = ci_cm.z;

    qlRe = (double*) malloc (((int)params[0] + 1) * sizeof (double));
    qlIm = (double*) malloc (((int)params[0] + 1) * sizeof (double));

    qlm((int)params[0], nk, xc, yc, zc, x, y, z, qlRe, qlIm);
    qss = qsum((int)params[0], qlRe, qlIm);

    double wlval = 0.0;

    if (qss > 1e-4)
    {
        wlval = Wl ((int)params[0], qlRe, qlIm) / (qss * qss * qss);
    }

    free (qlRe);
    free (qlIm);
    free (x);
    free (y);
    free (z);
    return wlval;
}

DerivedRegister<PackingWl> PackingWl::reg("PackingWl");