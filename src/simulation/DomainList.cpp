#include "DomainList.h"

DomainList::DomainList() : m(1), N(1), pbc(false), m_assigned(false), init_domains(false),
        x_min(0), y_min(0), z_min(0), x_max(0), y_max(0), z_max(0),
        dx(0), dy(0), dz(0), rc_max(0)
{}

DomainList::DomainList(const DomainList& orig) : m(orig.m), N(orig.N), pbc(orig.pbc), init_domains(orig.init_domains),
        m_assigned(orig.m_assigned), x_min(0), y_min(0), z_min(0), x_max(0), y_max(0), z_max(0), dx(0), dy(0), dz(0), rc_max(0) {}

DomainList::~DomainList() {}

void DomainList::setupDomainsList(double rcMax, Box& box)
{
    rc_max = rcMax;
    setBoxDim(box);
    initDomains();
    pbc = box.pbc;
}

void DomainList::initDomains()
{
    if (!init_domains)
    {
        int index;
        int neighix;
        for (int k = 0; k < m; k++)
        {
            for (int j = 0; j < m; j++)
            {
                for (int i = 0; i < m; i++)
                {
                    // ASSIGN INDEX
                    index = getDomainIndex(i, j, k);
                    domains[index].myid = index;
                
                    // EVERY DOMAIN IS ALSO ITS OWN NEIGHBOR
                    for (int l = -1; l <= 1; l++)
                      for (int o = -1; o <= 1; o++)
                        for (int p = -1; p <= 1; p++)
                        {
                            //neighix = getDomainIndex(i+v[l], j+v[o], k+v[p]);
                            neighix = getDomainIndex(i + l, j + o, k + p);
                            if (neighix != -1)
                            {
                                domains[index].addNeighDomain(neighix);
                            }
                      }
                }
            }
        }
        init_domains = true;
    }
    
}

int DomainList::getDomainIndex(int i, int j, int k)
{
    if (i < 0 && !pbc) {return -1;}
    else if (i  >= m  && !pbc) {return -1;}
    
    if (j < 0 && !pbc) {return -1;}
    else if (j  >= m  && !pbc) {return -1;}
    
    if (k < 0 && !pbc) {return -1;}
    else if (k  >= m  && !pbc) {return -1;}
    
    if (i < 0 && pbc) {i += m;}
    else if (i  >= m  && pbc) {i -= m;}
    
    if (j < 0 && pbc) {j += m;}
    else if (j  >= m  && pbc) {j -= m;}
    
    if (k < 0 && pbc) {k += m;}
    else if (k  >= m  && pbc) {k -= m;}
    
    return (i + j * m + k * m * m);
}

//void DomainList::assignCell(Cell& cell)
//{
//    
//    int cellid = cell.cellId;
//    //std::cout << "Assiging cell id=" << cellid << std::endl;
//    for (int i = 0; i < cell.numberofVertices(); i++)
//    {
//        assignVertex(cell.vertices[i], cellid);
//    }
//    //std::cout << "Assiging done" << std::endl;
//}

void DomainList::assignVertex(Vertex& vertex, int cellid)
{
    int index = getDomainIndex(vertex);
    domains[index].addVertex(vertex.getId(), cellid);
    vertex.domainIdx = index;
}

int DomainList::getDomainIndex(Vertex& vertex)
{
    int xix, yix, zix;
    
    double delx = vertex.xyz.x - x_min;
    double dely = vertex.xyz.y - y_min;
    double delz = vertex.xyz.z - z_min;
    
    xix = floor(delx / dx);
    yix = floor(dely / dy);
    zix = floor(delz / dz);
    return getDomainIndex(xix, yix, zix);
}

void DomainList::setBoxDim(Box& box)
{
    double lx, ly, lz;
    pbc = box.pbc;
    if (pbc)
    {    
        x_min = -box.getX();
        y_min = -box.getY();
        z_min = -box.getZ();
        
        x_max = box.getX();
        y_max = box.getY();
        z_max = box.getZ();
        
        //lx = 2 * box.getX();
        //ly = 2 * box.getY();
        //lz = 2 * box.getZ();
    } 
    else
    {
        x_min = -(box.getX() + rc_max);
        y_min = -(box.getY() + rc_max);
        z_min = -(box.getZ() + rc_max);
        
        x_max = box.getX() + rc_max;
        y_max = box.getY() + rc_max;
        z_max = box.getZ() + rc_max;
    }

    lx = x_max - x_min;
    ly = y_max - y_min;
    lz = z_max - z_min;

    
    if (!m_assigned)
    {
        double Lmax;
        Lmax = std::max(lx, ly);
        Lmax = std::max(lz, Lmax);
        m = ceil( Lmax / (2 * rc_max) );
        m = std::min(m, MAX_M);
        N = m * m * m;
        m_assigned = true;
    }
    
    dx = lx / m;
    dy = ly / m;
    dz = lz / m;
    
    //std::cout << "M=" << m << " number of domains=" << N << std::endl;
}

void DomainList::voidDomains()
{
    //std::cout << "czyscimy dupe" << std::endl;
    for (int i = 0; i < N; i++)
    {
        domains[i].voidParticlesInDomain();
    }
}
int DomainList::getNumberOfNeigh(int dix)
{
    return domains[dix].numofNeighDom;
}

int DomainList::getDomainNeighbor(int dix, int nix)
{
    return domains[dix].neighDomains[nix];
}

int DomainList::getVertexIdx(int dix, int xpart)
{
    return domains[dix].vertIds[xpart];
}

int DomainList::getCellIdx(int dix, int xpart)
{
    return domains[dix].cellsIds[xpart];
}

int DomainList::getNumOfParticles(int dix)
{
    return domains[dix].numofVert;
}

int DomainList::numberofAssignedParticles()
{
    int sum = 0;
    for (int i = 0; i < N; i++)
    {
        sum += domains[i].numofVert;
    }
    return sum;
}