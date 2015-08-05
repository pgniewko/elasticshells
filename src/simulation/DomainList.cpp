#include "DomainList.h"

utils::Logger DomainList::domainlist_logs("domainlist");

DomainList::DomainList() : m(1), N(1), pbc(false), m_assigned(false), init_domains(false),
    x_min(0), y_min(0), z_min(0),
    x_max(0), y_max(0), z_max(0),
    dx(0), dy(0), dz(0), rc_max(0) 
{
}

DomainList::DomainList(const DomainList& orig) : m(orig.m), N(orig.N), pbc(orig.pbc),
    m_assigned(orig.m_assigned), init_domains(orig.init_domains),
    x_min(0), y_min(0), z_min(0),
    x_max(0), y_max(0), z_max(0),
    dx(0), dy(0), dz(0), rc_max(0) 
{
    for (int i = 0; i < N; i++)
        vertsInDomains[i] = orig.vertsInDomains[i];
    
    std::cout << "SAFETY CHECK" << std::endl;
    exit(EXIT_FAILURE);
}

DomainList::~DomainList() 
{
}

void DomainList::setupDomainsList(double rcMax, Box& box)
{
    rc_max = rcMax;
    pbc = box.pbc;

    if (!m_assigned)
    {
        setM(box);
    }

    setBoxDim(box);
    initDomains();
}

void DomainList::initDomains()
{
    if (!init_domains)
    {
        int index;
        int neighix = -1;

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
                                neighix = getDomainIndex(i + l, j + o, k + p);

                                if (neighix != -1)
                                {
                                    addNeighDomain(index, neighix);
                                    //domains[index].addNeighDomain(neighix);
                                }
                            }
                }
            }
        }

        init_domains = true;
        domainlist_logs << utils::LogLevel::FINEST << "Linked domains list has been initialized successfully." << "\n";
    }
}

void DomainList::addNeighDomain(int dix, int nidx)
{
    try
    {
        //std::cout << domains[dix].numberOfNeighs << " " << MAX_D_NEIGH << std::endl;
        if (domains[dix].numberOfNeighs >= MAX_D_NEIGH)
            throw MaxSizeException("Trying to add more domain neighbors than it's possible.\n"
                                   "New neighbor will not be added!\n"
                                   "The code will be terminated due to the bug!\n");

        domains[dix].neighborDomainIdx[ domains[dix].numberOfNeighs ] = nidx;
        domains[dix].numberOfNeighs++;
      
        //neighborDomainIdx[numberOfNeighs] = domId;
        //numberOfNeighs++;
    }
    catch (MaxSizeException& e)
    {
        domainlist_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        exit(EXIT_FAILURE);
    }    
//    int num_nieg = domains[dix].numberOfNeighs;
//    domains[dix].neighborDomainIdx[ domains[dix].numberOfNeighs ] = nidx;
//    domains[dix].numberOfNeighs++;
}

int DomainList::getDomainIndex(int i, int j, int k)
{
    if (i < 0 && !pbc)
    {
        return -1;
    }
    else if (i  >= m  && !pbc)
    {
        return -1;
    }

    if (j < 0 && !pbc)
    {
        return -1;
    }
    else if (j  >= m  && !pbc)
    {
        return -1;
    }

    if (k < 0 && !pbc)
    {
        return -1;
    }
    else if (k  >= m  && !pbc)
    {
        return -1;
    }

    if (i < 0 && pbc)
    {
        i += m;
    }
    else if (i  >= m  && pbc)
    {
        i -= m;
    }

    if (j < 0 && pbc)
    {
        j += m;
    }
    else if (j  >= m  && pbc)
    {
        j -= m;
    }

    if (k < 0 && pbc)
    {
        k += m;
    }
    else if (k  >= m  && pbc)
    {
        k -= m;
    }

    return (i + j * m + k * m * m);
}

void DomainList::assignVertex(Vertex& vertex, int cellid)
{
    int index = getDomainIndex(vertex);
//    domains[index].addVertex(vertex.getId(), cellid);
//    vertex.domainIdx = index;
    
    try
    {
        if (vertsInDomains[index] >= MAX_IN_DOMAIN)
            throw MaxSizeException("Trying to add more vertices than it's allowed.\n"
                                   "This may significantly affect the simulation accuracy!\n"
                                   "Simulation is about to end.");

        
        domains[index].vertIds[ vertsInDomains[index] ] = vertex.getId();
        domains[index].cellsIds[ vertsInDomains[index] ] = cellid;
        vertsInDomains[index]++;
        vertex.domainIdx = index;
        //std::cout << "Assiging a vertex="<<vertex.getId() << " cellid="<<cellid << " vertsInDomains["<<index<<"]="<<vertsInDomains[index]<< std::endl;
    }
    catch (MaxSizeException& e)
    {
        domainlist_logs << utils::LogLevel::INFO << "MY_ID=" << domains[index].myid << "\n";
        domainlist_logs << utils::LogLevel::INFO << "NUMBER_OF_VERTICES=" << vertsInDomains[index] << "\n";
        domainlist_logs << utils::LogLevel::INFO << "TRYING_TO_ADD:VERTEX_ID=" << vertex.getId() << " CELL_ID=" << cellid << "\n";
        domainlist_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
    
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

void DomainList::setM(Box& box)
{
    double lx, ly, lz;
    double xmin, xmax, ymin, ymax, zmin, zmax;

    if (pbc)
    {
        // MAKE SURE THAT AT THE END of SIM. DOMAINS ARE NOT TOO SMALL !!!
        xmin = -box.getXend();
        ymin = -box.getYend();
        zmin = -box.getZend();
        xmax =  box.getXend();
        ymax =  box.getYend();
        zmax =  box.getZend();
    }
    else
    {
        xmin = -(box.getXend() + rc_max);
        ymin = -(box.getYend() + rc_max);
        zmin = -(box.getZend() + rc_max);
        xmax =   box.getXend() + rc_max;
        ymax =   box.getYend() + rc_max;
        zmax =   box.getZend() + rc_max;
    }

    lx = xmax - xmin;
    ly = ymax - ymin;
    lz = zmax - zmin;
    double Lmax;
    Lmax = std::max(lx, ly);
    Lmax = std::max(lz, Lmax);
    m = ceil( Lmax / (2 * rc_max + EPSILON) );
    m = std::min(m, MAX_M);
    N = m * m * m;
    m_assigned = true;
    domains.reserve(N);
    domainlist_logs << utils::LogLevel::INFO << "NUMBER OF LINKED DOMAINS N_DOMAINS=" << N << "\n";
}

void DomainList::setBoxDim(Box& box)
{
    double lx, ly, lz;

    if (pbc)
    {
        x_min = -box.getX();
        y_min = -box.getY();
        z_min = -box.getZ();
        x_max =  box.getX();
        y_max =  box.getY();
        z_max =  box.getZ();
    }
    else
    {
        x_min = -(box.getX() + rc_max);
        y_min = -(box.getY() + rc_max);
        z_min = -(box.getZ() + rc_max);
        x_max =   box.getX() + rc_max;
        y_max =   box.getY() + rc_max;
        z_max =   box.getZ() + rc_max;
    }

    lx = x_max - x_min;
    ly = y_max - y_min;
    lz = z_max - z_min;
    dx = lx / m;
    dy = ly / m;
    dz = lz / m;
}

void DomainList::voidDomains()
{
    for (int i = 0; i < N; i++)
        vertsInDomains[i] = 0;
    
    //for (int i = 0; i < N; i++)
    //{
    //    domains[i].voidParticlesInDomain();
    //}
}
int DomainList::getNumberOfNeigh(int dix)
{
    return domains[dix].numberOfNeighs;
}

int DomainList::getDomainNeighbor(int dix, int nix)
{
    return domains[dix].neighborDomainIdx[nix];
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
    //return domains[dix].numberOfVerts;
    return vertsInDomains[dix];
}

int DomainList::numberofAssignedParticles()
{
    int sum = 0;

    for (int i = 0; i < N; i++)
        sum += vertsInDomains[i]; 
    
//    for (int i = 0; i < N; i++)
//    {
//        sum += domains[i].numberOfVerts;
//    }

    return sum;
}

double DomainList::getMaxScale()
{
    return rc_max;
}