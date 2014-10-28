#include "Domain.h"

utils::Logger domain_logs("domain");

Domain::Domain() : myid(-1), numberOfNeighs(0), numberOfVerts(0) {}

Domain::Domain(const Domain& orig) : myid(orig.myid),  numberOfNeighs(orig.numberOfNeighs), numberOfVerts(orig.numberOfVerts)
{
    for (int i = 0; i < numberOfVerts; i++)
    {
        vertIds[i] = orig.vertIds[i];
        cellsIds[i] = orig.cellsIds[i];
    }
    for (int i = 0; i < numberOfNeighs; i++)
    {
        neighborDomainIdx[i] = orig.neighborDomainIdx[i];
    }
}

Domain::~Domain() {}

void Domain::addVertex(int vid, int cellid)
{
    try {
        if (numberOfVerts >= MAX_IN_DOMAIN)
            throw MaxSizeException("Trying to add more vertices than it's allowed."
                                   "New vertex will not be added !"
                                   "This may significantly affect the simulation accuracy !");
        
        vertIds[numberOfVerts] = vid;
        cellsIds[numberOfVerts] = cellid;
        numberOfVerts++;
    }
    catch (MaxSizeException& e)
    {
        domain_logs << utils::LogLevel::WARNING << e.what() << "\n";
        return;
    }
}

void Domain::addNeighDomain(int domId)
{
    try {
        if (numberOfVerts >= MAX_IN_DOMAIN)
            throw MaxSizeException("Trying to add more domain neighbors than it's possible."
                                   "New neighbor will not be added!"
                                   "The code will be terminated due the bug !");
        
        neighborDomainIdx[numberOfNeighs] = domId;
        numberOfNeighs++;
    } 
    catch (MaxSizeException& e)
    {
        domain_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        exit(1);
    }
}
void Domain::voidParticlesInDomain()
{
    numberOfVerts = 0;
}