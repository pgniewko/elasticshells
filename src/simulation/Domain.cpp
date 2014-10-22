#include "Domain.h"

Domain::Domain() : myid(-1), numofVert(0), numofNeighDom(0) {}

Domain::Domain(const Domain& orig) : myid(orig.myid), numofVert(orig.numofVert), numofNeighDom(orig.numofNeighDom) 
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
    if (numofVert < MAX_IN_DOMAIN)
    {
        vertIds[numofVert] = vid;
        cellsIds[numofVert] = cellid;
        numofVert++;
    }
    else
    {
        std::cout << "dupa jasia" << std::endl;
        exit(1);
    } 
}

void Domain::addNeighDomain(int domId)
{
    if (numofNeighDom < 27)
    {
        neighDomains[numofNeighDom] = domId;
        numofNeighDom++;
    }
    else
    {
        std::cout << "dupa zosi" << std::endl;
        exit(1);
        // exception
    }
}
void Domain::voidParticlesInDomain()
{
    numofVert = 0;
}


