#include "SurfacePressure.h"

utils::Logger SurfacePressure::sp_log("surface_press");

SurfacePressure::SurfacePressure(const char* name, const char* format) : Observer(name, format) {}

SurfacePressure::SurfacePressure(const SurfacePressure& orig) : Observer(orig) {}

SurfacePressure::~SurfacePressure() {}

void SurfacePressure::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
    d_param = strtod(args_[ num + 1 ].c_str(), NULL);
}

double SurfacePressure::observe(const Box& box, std::vector<Cell>& cells)
{
    double pressure = 0.0;
    
    if ( i_param == 0 )
    {
        double r1;
        double r2;
        double e1;
        double e2;
        double nu1;
        double nu2;
              
        double volume = box.getVolume(d_param);
        Vector3D rij;
        Vector3D fij;
        
        std::size_t n = cells.size();
        for (std::size_t k = 0; k < n; k++)
        {
            for (std::size_t l = k+1; l < n; l++)
            {
                r1 = cells[k].getVertexR();
                r2 = cells[l].getVertexR();
                e1 = cells[k].getE();
                e2 = cells[l].getE();
                nu1 = cells[k].getNu();
                nu2 = cells[l].getNu();
                
                for (int i = 0; i < cells[k].getNumberVertices(); i++ )
                {
                    for (int j = 0; j < cells[l].getNumberVertices(); j++ )
                    {
                         Box::getDistance(rij, cells[k].vertices[i].r_c, cells[l].vertices[j].r_c, box);
                         fij = HertzianRepulsion::calcForce(rij, r1, r2, e1, e2, nu1, nu2);
                         pressure += dot(rij, fij);
                    }
                }
               
            }
        }
        
        // Two other terms missing:
        // i) box-vertex interaction
        // ii) strains in the finite-elements
        
        pressure /= (3.0*volume);
        return pressure;
    }
    
    else if ( i_param == 1)
    {
        if (box.pbc)
        {
            sp_log << utils::LogLevel::INFO << "\n";
            return 0.0;
        }
        
        double totalForce = SurfaceForce::calcTotalForce(box, cells);
        double area = box.getArea(d_param);
        pressure = totalForce / area;    
    }
    else
    {
        return 0.0;
    }
    
    return pressure;
}

DerivedRegister<SurfacePressure> SurfacePressure::reg("SurfacePressure");