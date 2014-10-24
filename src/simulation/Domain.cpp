#include "Domain.h"

utils::Logger domain_logs("domain");

Domain::Domain() : myid(-1), numofNeighDom(0), numofVert(0) {}

Domain::Domain(const Domain& orig) : myid(orig.myid),  numofNeighDom(orig.numofNeighDom), numofVert(orig.numofVert)
{
    for (int i = 0; i < numofVert; i++)
    {
        vertIds[i] = orig.vertIds[i];
        cellsIds[i] = orig.cellsIds[i];
    }
    for (int i = 0; i < numofNeighDom; i++)
    {
        neighDomains[i] = orig.neighDomains[i];
    }
}

Domain::~Domain() {}

void Domain::addVertex(int vid, int cellid)
{
    try {
        
    
        if (numofVert >= MAX_IN_DOMAIN)
            throw MaxSizeException("Trying to add more vertices than it's allowed."
                                   "New vertex will not be added !"
                                   "This may significantly affect the simulation accuracy !");
        
        vertIds[numofVert] = vid;
        cellsIds[numofVert] = cellid;
        numofVert++;
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
        if (numofVert >= MAX_IN_DOMAIN)
            throw MaxSizeException("Trying to add more domain neighbors than it's possible."
                                   "New neighbor will not be added!"
                                   "The code will be terminated due the bug !");
        
        neighDomains[numofNeighDom] = domId;
        numofNeighDom++;
    } 
    catch (MaxSizeException& e)
    {
        domain_logs << utils::LogLevel::CRITICAL << e.what() << "\n";
        exit(1);
    }
}
void Domain::voidParticlesInDomain()
{
    numofVert = 0;
}


