#ifndef DOMAIN_H
#define	DOMAIN_H

#include "Environment.h"

class Domain {
public:
    Domain();
    Domain(const Domain& orig);
    virtual ~Domain();
    
    void addVertex(int, int);
    void addNeighDomain(int);
    void voidParticlesInDomain();
    int vertIds[MAX_IN_DOMAIN];
    int cellsIds[MAX_IN_DOMAIN];
    int neighDomains[27];
    
    int myid;
    int numofNeighDom;
    int numofVert;
    
private:

    //double lx, ly, lz;
    //double x0, y0, z0;
};

#endif	/* DOMAIN_H */

