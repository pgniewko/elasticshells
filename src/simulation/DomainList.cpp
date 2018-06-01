#include "DomainList.h"

utils::Logger DomainList::domainlist_logs("domainlist");

DomainList::DomainList() : m(1), N(1), pbc(false), m_assigned(false), init_domains(false),
    x_min(0), y_min(0), z_min(0),
    x_max(0), y_max(0), z_max(0),
    dx(0), dy(0), dz(0), minL(0)
{
}

DomainList::DomainList(const DomainList& orig) : m(orig.m), N(orig.N), pbc(orig.pbc),
    m_assigned(orig.m_assigned), init_domains(orig.init_domains),
    x_min(0), y_min(0), z_min(0),
    x_max(0), y_max(0), z_max(0),
    dx(0), dy(0), dz(0), minL(0)
{
    for (int i = 0; i < N; i++)
    {
        head[i] = NULL;
    }

    //TODO: implement proper copying.
    domainlist_logs <<  utils::LogLevel::SEVERE << "THIS CLASS SHOULD NEVER BE COPIED!\n";
    domainlist_logs <<  utils::LogLevel::SEVERE << "PROGRAM TERMINATION FOR SAFETY REASONS!\n";
    exit(EXIT_FAILURE);
}

DomainList::~DomainList() {}

void DomainList::setupDomainsList(double lengthScale, Box& box)
{
    minL = lengthScale + constants::delta7;
    pbc = box.pbc;

    if (!m_assigned)
    {
        setM(box);
    }

    initDomains();
    setBoxDim(box);
}

bool DomainList::validateLinkedDomains()
{
    for (int i = 0; i < N; i++)
    {
        if (domains[i].myid < 0)
        {
            return false;
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < domains[i].neighborDomainNumber; j++)
        {
            if (domains[i].neighborDomainIdx[j] < 0)
            {
                return false;
            }
        }
    }

    return true;
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

                    for (int l = -1; l <= 1; l++)
                        for (int o = -1; o <= 1; o++)
                            for (int p = -1; p <= 1; p++)
                            {
                                neighix = getDomainIndex(i + l, j + o, k + p);

                                if (neighix != -1 && neighix != index)
                                {
                                    addNeighDomain(index, neighix);
                                }
                            }
                }
            }
        }

        if ( ! validateLinkedDomains() )
        {
            domainlist_logs << utils::LogLevel::SEVERE << "DOMAIN INITIALIZATION HAS FAILED!\n";
            exit(EXIT_FAILURE);
        }
        else
        {
            init_domains = true;
            domainlist_logs << utils::LogLevel::FINEST << "Linked domains list has been initialized successfully." << "\n";
        }
    }
}

void DomainList::addNeighDomain(int dix, int nidx)
{
    try
    {
        if (domains[dix].neighborDomainNumber >= MAX_D_NEIGH)
        {
            throw MaxSizeException("Trying to add more domain neighbors than it's possible.\n"
                                   "New neighbor will not be added!\n"
                                   "The code will be terminated due to the BUG!\n");
        }

        domains[dix].neighborDomainIdx[ domains[dix].neighborDomainNumber ] = nidx;
        domains[dix].neighborDomainNumber++;
    }
    catch (MaxSizeException& e)
    {
        domainlist_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        domainlist_logs << "domains[" << dix << "].neighborDomainNumber=" << domains[dix].neighborDomainNumber << "\n";
        exit(EXIT_FAILURE);
    }
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

void DomainList::assignVertex(Vertex* vertex)
{
    int index = getDomainIndex(vertex);
    vertex->next = head[index];
    head[index] = vertex;
}

int DomainList::getDomainIndex(Vertex* vertex)
{
    int xix, yix, zix;
    double delx = vertex->r_c.x - x_min;
    double dely = vertex->r_c.y - y_min;
    double delz = vertex->r_c.z - z_min;
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
        xmin = -box.getX();
        ymin = -box.getY();
        zmin = -box.getZ();
        xmax =  box.getX();
        ymax =  box.getY();
        zmax =  box.getZ();
    }
    else
    {
        xmin = -(box.getX() + 0.5*minL);
        ymin = -(box.getY() + 0.5*minL);
        zmin = -(box.getZ() + 0.5*minL);
        xmax =   box.getX() + 0.5*minL;
        ymax =   box.getY() + 0.5*minL;
        zmax =   box.getZ() + 0.5*minL;
    }

    lx = xmax - xmin;
    ly = ymax - ymin;
    lz = zmax - zmin;
    double Lmax = std::max( std::max(lx, ly), lz );
    m = ceil( Lmax / minL );
    m = std::min(m, MAX_M);
    N = m * m * m;
    m_assigned = true;
    domains.reserve(N);

    for (int i = 0; i < N; i++)
    {
        domains[i] = domain_t();
    }

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
        x_min = -(box.getX() + 0.5*minL);
        y_min = -(box.getY() + 0.5*minL);
        z_min = -(box.getZ() + 0.5*minL);
        x_max =   box.getX() + 0.5*minL;
        y_max =   box.getY() + 0.5*minL;
        z_max =   box.getZ() + 0.5*minL;
    }

    lx = x_max - x_min;
    ly = y_max - y_min;
    lz = z_max - z_min;
    dx = lx / m;
    dy = ly / m;
    dz = lz / m;
    
    dx = std::min(dx, minL);
    dy = std::min(dy, minL);
    dz = std::min(dz, minL);    
}

void DomainList::voidDomains()
{
    for (int i = 0; i < N; i++)
    {
        head[i] = 0;
    }
}
int DomainList::getNumberOfNeigh(int dix)
{
    return domains[dix].neighborDomainNumber;
}

int DomainList::getNumOfParticles(int dix)
{
    int counter = 0;

    for (Vertex* s = head[dix]; s != 0; s = s->next)
    {
        counter++;
    }

    return counter;
}

int DomainList::numberofAssignedParticles()
{
    int sum = 0;

    for (int i = 0; i < N; i++)
    {
        sum += getNumOfParticles(i);
    }

    return sum;
}

double DomainList::getMinLength()
{
    return minL;
}

void DomainList::calcNbForces(std::vector<Shell>& cells, const Box& box) const
{
    Vertex* target;
    Vertex* partner;

    int neighIndex;

    for (int domainIdx = 0; domainIdx < N; domainIdx++)
    {
        if (head[domainIdx] != 0)
        {
            // INTRA-DOMAIN CONTACTS
            for (target = head[domainIdx]; target != 0; target = target->next)
            {
                for (partner = target->next; partner != 0; partner = partner->next)
                {
                    nbForce(target, partner, cells, box);
                }
            }

            // INTER-DOMAIN CONTACTS
            for (target = head[domainIdx]; target != 0; target = target->next)
            {
                for (int k = 0; k < domains[domainIdx].neighborDomainNumber; k++)
                {
                    neighIndex = domains[domainIdx].neighborDomainIdx[k];

                    if (neighIndex > domainIdx)
                    {
                        for (partner = head[neighIndex]; partner != 0; partner = partner->next)
                        {
                            nbForce(target, partner, cells, box);
                        }
                    }
                }
            }
        }
    }
}


void DomainList::nbForce(Vertex* target, Vertex* partner, std::vector<Shell>& cells, const Box& box) const
{
    int cellId_target = target->getCellId();
    int cellId_partner = partner->getCellId();

    const struct shell_params_t params1 = cells[cellId_target].get_params();
    const struct shell_params_t params2 = cells[cellId_partner].get_params();

    double r1 = params1.vertex_r;
    double r2 = params2.vertex_r;
    double e1 = params1.ecc;
    double e2 = params2.ecc;
    double nu1 = params1.nu;
    double nu2 = params2.nu;

    Vector3D force(0, 0, 0);
    Vector3D dij;

    if (cellId_target != cellId_partner)
    {
        Box::getDistance(dij, partner->r_c, target->r_c, box);
        force = HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
        target->f_c  += force;
        partner->f_c += -force;
    }
    else
    {
        int i = target->getId();
        int j = partner->getId();

        if (i != j && !target->isNeighbor(j))
        {
            Box::getDistance(dij, partner->r_c, target->r_c, box);
            force = HertzianRepulsion::calcForce(dij, r1, r1, e1, e1, nu1, nu1);
            target->f_c  += force;
            partner->f_c += -force;
        }
    }
}


// *****************************************************************************
Vector3D DomainList::getNbForce(Vertex* target, Vertex* partner, const std::vector<Shell>& cells, const Box& box) const
{

    int cellId_target = target->getCellId();
    int cellId_partner = partner->getCellId();

    const struct shell_params_t params1 = cells[cellId_target].get_params();
    const struct shell_params_t params2 = cells[cellId_partner].get_params();

    double r1 = params1.vertex_r;
    double r2 = params2.vertex_r;
    double e1 = params1.ecc;
    double e2 = params2.ecc;
    double nu1 = params1.nu;
    double nu2 = params2.nu;

    Vector3D force(0, 0, 0);
    Vector3D dij;

    if (cellId_target != cellId_partner)
    {
        Box::getDistance(dij, partner->r_c, target->r_c, box);
        force = HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);
    }

//    else
//    {
//        int i = target->getId();
//        int j = partner->getId();
//
//        if (i != j && !target->isNeighbor(j))
//        {
//            Box::getDistance(dij, partner->r_c, target->r_c, box);
//            force = HertzianRepulsion::calcForce(dij, r1, r1, e1, e1, nu1, nu1);
//        }
//    }

    return force;
}

double DomainList::virial(Vertex* target, Vertex* partner, const std::vector<Shell>& cells, const Box& box) const
{
    int cellId_target = target->getCellId();
    int cellId_partner = partner->getCellId();

    const struct shell_params_t params1 = cells[cellId_target].get_params();
    const struct shell_params_t params2 = cells[cellId_partner].get_params();

    double r1 = params1.vertex_r;
    double r2 = params2.vertex_r;
    double e1 = params1.ecc;
    double e2 = params2.ecc;
    double nu1 = params1.nu;
    double nu2 = params2.nu;

    Vector3D fij(0, 0, 0);
    Vector3D dij;

    double pij = 0.0;

    if (cellId_target != cellId_partner)
    {
        Box::getDistance(dij, partner->r_c, target->r_c, box);
        fij = HertzianRepulsion::calcForce(dij, r1, r2, e1, e2, nu1, nu2);

        if (fij.length_sq() > Shell::MIN_FORCE_SQ)
        {
            pij = dot(dij, fij);
        }
    }
    else
    {
        int i = target->getId();
        int j = partner->getId();

        if (i != j && !target->isNeighbor(j))
        {
            Box::getDistance(dij, partner->r_c, target->r_c, box);
            fij = HertzianRepulsion::calcForce(dij, r1, r1, e1, e1, nu1, nu1);

            if (fij.length_sq() > Shell::MIN_FORCE_SQ)
            {
                pij = dot(dij, fij);
            }
        }
    }

    return pij;
}

double DomainList::calcContactForce(const int cell1id, const int cell2id, const std::vector<Shell>& cells, const Box& box) const
{
    Vertex* target;
    Vertex* partner;

    int neighIndex;

    Vector3D force;

    std::vector<Vector3D> verts_forces;

    for (int vix = 0; vix < cells[cell1id].getNumberVertices(); vix++)
    {
        Vector3D new_one(0, 0, 0);
        verts_forces.push_back(  new_one );
    }


    double contact_force = 0.0;

    for (int domainIdx = 0; domainIdx < N; domainIdx++)
    {
        if (head[domainIdx] != 0)
        {
            // INTRA-DOMAIN CONTACTS
            for (target = head[domainIdx]; target != 0; target = target->next)
            {
                for (partner = target->next; partner != 0; partner = partner->next)
                {
                    if ( target->getCellId() == cell1id && partner->getCellId() == cell2id )
                    {
                        force = getNbForce(target, partner, cells, box);
                        verts_forces[ target->getId() ] += force;
                    }
                    else if ( target->getCellId() == cell2id && partner->getCellId() == cell1id )
                    {
                        force = getNbForce(partner, target, cells, box);
                        verts_forces[ partner->getId() ] += force;
                    }
                }
            }

            // INTER-DOMAIN CONTACTS
            for (target = head[domainIdx]; target != 0; target = target->next)
            {
                for (int k = 0; k < domains[domainIdx].neighborDomainNumber; k++)
                {
                    neighIndex = domains[domainIdx].neighborDomainIdx[k];

                    if (neighIndex > domainIdx)
                    {
                        for (partner = head[neighIndex]; partner != 0; partner = partner->next)
                        {
                            if ( target->getCellId() == cell1id && partner->getCellId() == cell2id )
                            {
                                force = getNbForce(target, partner, cells, box);
                                verts_forces[ target->getId() ] += force;
                            }
                            else if ( target->getCellId() == cell2id && partner->getCellId() == cell1id )
                            {
                                force = getNbForce(partner, target, cells, box);
                                verts_forces[ partner->getId() ] += force;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int vix = 0; vix < cells[cell1id].getNumberVertices(); vix++)
    {
        if (verts_forces[vix].length_sq() > 0.0 )
        {
            contact_force += verts_forces[vix].length();
        }
    }

    return contact_force;
}

double DomainList::virialPressure(const int cell1id, const int cell2id, const std::vector<Shell>& cells, const Box& box) const
{
    double pressure = 0.0;

    Vertex* target;
    Vertex* partner;

    int neighIndex;

    for (int domainIdx = 0; domainIdx < N; domainIdx++)
    {
        if (head[domainIdx] != 0)
        {
            // INTRA-DOMAIN CONTACTS
            for (target = head[domainIdx]; target != 0; target = target->next)
            {
                for (partner = target->next; partner != 0; partner = partner->next)
                {
                    if ( (target->getCellId() == cell1id && partner->getCellId() == cell2id) || (target->getCellId() == cell2id && partner->getCellId() == cell1id) )
                    {
                        pressure += virial(target, partner, cells, box);
                    }
                }
            }

            // INTRA-DOMAIN CONTACTS
            for (target = head[domainIdx]; target != 0; target = target->next)
            {
                for (int k = 0; k < domains[domainIdx].neighborDomainNumber; k++)
                {
                    neighIndex = domains[domainIdx].neighborDomainIdx[k];

                    if (neighIndex > domainIdx)
                    {
                        for (partner = head[neighIndex]; partner != 0; partner = partner->next)
                        {
                            if ( (target->getCellId() == cell1id && partner->getCellId() == cell2id) || (target->getCellId() == cell2id && partner->getCellId() == cell1id) )
                            {
                                pressure += virial(target, partner, cells, box);
                            }
                        }
                    }
                }
            }
        }
    }

    return pressure;
}


bool DomainList::isInContact(const int cell1id, const int cell2id, const std::vector<Shell>& cells, const Box& box) const
{
    Vertex* target;
    Vertex* partner;

    int neighIndex;

    Vector3D force;


    for (int domainIdx = 0; domainIdx < N; domainIdx++)
    {
        if (head[domainIdx] != 0)
        {
            // INTRA-DOMAIN CONTACTS
            for (target = head[domainIdx]; target != 0; target = target->next)
            {
                for (partner = target->next; partner != 0; partner = partner->next)
                {
                    if ( (target->getCellId() == cell1id && partner->getCellId() == cell2id) || (target->getCellId() == cell2id && partner->getCellId() == cell1id) )
                    {
                        force = getNbForce(target, partner, cells, box);

                        if (force.length_sq() > Shell::MIN_FORCE_SQ)
                        {
                            return true;
                        }
                    }
                }
            }

            // INTRA-DOMAIN CONTACTS
            for (target = head[domainIdx]; target != 0; target = target->next)
            {
                for (int k = 0; k < domains[domainIdx].neighborDomainNumber; k++)
                {
                    neighIndex = domains[domainIdx].neighborDomainIdx[k];

                    for (partner = head[neighIndex]; partner != 0; partner = partner->next)
                    {
                        if ( (target->getCellId() == cell1id && partner->getCellId() == cell2id) || (target->getCellId() == cell2id && partner->getCellId() == cell1id) )
                        {
                            force = getNbForce(target, partner, cells, box);

                            if (force.length_sq() > Shell::MIN_FORCE_SQ)
                            {
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }

    return false;
}

double DomainList::calcNbEnergy(const std::vector<Shell>& cells, const Box& box) const
{
    double totalNbEnergy = 0.0;
    Vertex* target;
    Vertex* partner;

    int neighIndex;

    for (int domainIdx = 0; domainIdx < N; domainIdx++)
    {
        if (head[domainIdx] != 0)
        {
            // INTRA-DOMAIN CONTACTS
            for (target = head[domainIdx]; target != 0; target = target->next)
            {
                for (partner = target->next; partner != 0; partner = partner->next)
                {
                    totalNbEnergy += nbEnergy(target, partner, cells, box);
                }
            }

            // INTRA-DOMAIN CONTACTS
            for (target = head[domainIdx]; target != 0; target = target->next)
            {
                for (int k = 0; k < domains[domainIdx].neighborDomainNumber; k++)
                {
                    neighIndex = domains[domainIdx].neighborDomainIdx[k];

                    if (neighIndex > domainIdx)
                    {
                        for (partner = head[neighIndex]; partner != 0; partner = partner->next)
                        {
                            totalNbEnergy += nbEnergy(target, partner, cells, box);
                        }
                    }
                }
            }
        }
    }

    return totalNbEnergy;
}

double DomainList::nbEnergy(const Vertex* target, const Vertex* partner, const std::vector<Shell>& cells, const Box& box) const
{
    double nb_energy = 0.0;
    int cellId_target = target->getCellId();
    int cellId_partner = partner->getCellId();

    const struct shell_params_t params1 = cells[cellId_target].get_params();
    const struct shell_params_t params2 = cells[cellId_partner].get_params();

    double r1 = params1.vertex_r;
    double r2 = params2.vertex_r;
    double e1 = params1.ecc;
    double e2 = params2.ecc;
    double nu1 = params1.nu;
    double nu2 = params2.nu;

    Vector3D dij;

    if (cellId_target != cellId_partner)
    {
        Box::getDistance(dij, partner->r_c, target->r_c, box);
        nb_energy = HertzianRepulsion::calcEnergy(dij, r1, r2, e1, e2, nu1, nu2);
    }
    else
    {
        int i = target->getId();
        int j = partner->getId();

        if (i != j && !target->isNeighbor(j))
        {
            Box::getDistance(dij, partner->r_c, target->r_c, box);
            nb_energy = HertzianRepulsion::calcEnergy(dij, r1, r1, e1, e1, nu1, nu1);

        }
    }

    return nb_energy;
}