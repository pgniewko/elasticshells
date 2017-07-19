#include "Roundness.h"

Roundness::Roundness(const char* name, const char* format) : Observer(name, format) {}

Roundness::Roundness(const Roundness& orig) : Observer(orig) {}

Roundness::~Roundness() {}

void Roundness::set_params(const int num, std::vector<std::string> args_)
{
    return;
};

double Roundness::observe(const Box& box, std::vector<Cell>& cells, const DomainList& dl)
{
    double V, D, DR, roundness_ = 0.0;
    
    for (uint i = 0; i < cells.size(); i++)
    {
        V = cells[i].calcVolume();
        V = 3.0 * V / (4.0 * constants::pi);
        D = 2.0 * std::pow( V, 1.0/3.0 ); // HEYWOOD'S DIAMETER
        DR = 2.0 * miniball_r(cells[i]);
        
        if (D > 0)
        {
            roundness_ += ( DR / D );
        }
        else
        {
            throw DataException("Division by zero");
        }
    }
    
    roundness_ /= cells.size();
            
    return roundness_;
}

double Roundness::miniball_r(Cell& cell)
{
    typedef double mytype;            // coordinate type

    int             d = 3;            // dimension
    int             n;      // number of points
    
    n = cell.getNumberVertices();
    
    mytype** ap = new mytype*[n];
    
    for (int i = 0; i < n; ++i)
    {
        mytype* p = new mytype[d];

        p[0] = cell.vertices[i].r_c.x;
        p[1] = cell.vertices[i].r_c.y;
        p[2] = cell.vertices[i].r_c.z;
        
        ap[i] = p;
    }
    
    // define the types of iterators through the points and their coordinates
    // ----------------------------------------------------------------------
    typedef mytype* const* PointIterator; 
    typedef const mytype* CoordIterator;

    // create an instance of Miniball
    // ------------------------------
    typedef Miniball::Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > MB;
  
    MB mb (d, ap, ap+n);
    
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