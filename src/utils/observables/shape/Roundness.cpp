#include "Roundness.h"

Roundness::Roundness(const char* name, const char* format) : Observer(name, format) {}

Roundness::Roundness(const Roundness& orig) : Observer(orig) {}

Roundness::~Roundness() {}

double Roundness::observe(const Box& box, const std::vector<Shell>& shells)
{
    double V, D, DR, roundness_ = 0.0;

    for (uint i = 0; i < shells.size(); i++)
    {
        V = shells[i].calc_volume();
        V = 3.0 * V / (4.0 * constants::pi);
        D = 2.0 * std::pow( V, 1.0 / 3.0 ); // HEYWOOD'S DIAMETER
        DR = 2.0 * miniball_r(shells[i]);

        if (D > 0)
        {
            roundness_ += ( DR / D );
        }
        else
        {
            throw DataException("Division by zero");
        }
    }

    roundness_ /= shells.size();

    return roundness_;
}

double Roundness::miniball_r(const Shell& shell)
{
    typedef double mytype;            // coordinate type

    int             d = 3;            // dimension
    int             n;      // number of points

    n = shell.get_number_vertices();

    mytype** ap = new mytype*[n];

    for (int i = 0; i < n; ++i)
    {
        mytype* p = new mytype[d];

        p[0] = shell.vertices[i].r_c.x;
        p[1] = shell.vertices[i].r_c.y;
        p[2] = shell.vertices[i].r_c.z;

        ap[i] = p;
    }

    // define the types of iterators through the points and their coordinates
    // ----------------------------------------------------------------------
    typedef mytype* const* PointIterator;
    typedef const mytype* CoordIterator;

    // create an instance of Miniball
    // ------------------------------
    typedef Miniball::Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > MB;

    MB mb (d, ap, ap + n);

    // squared radius
    mytype r2 =  mb.squared_radius();
    mytype r = std::sqrt( (double) r2 );

    // clean up
    for (int i = 0; i < n; ++i)
    {
        delete[] ap[i];
    }

    delete[] ap;

    return (double) r;
}

DerivedRegister<Roundness> Roundness::reg("Roundness");