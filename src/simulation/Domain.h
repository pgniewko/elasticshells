#ifndef DOMAIN_H
#define	DOMAIN_H


#include "Environment.h"
#include "exceptions/MaxSizeException.h"

class Domain
{
    public:
        Domain();
        Domain(const Domain& orig);
        virtual ~Domain();

        void addVertex(int, int);
        void addNeighDomain(int);
        void voidParticlesInDomain();
        int vertIds[MAX_IN_DOMAIN];
        int cellsIds[MAX_IN_DOMAIN];
        int neighborDomainIdx[27];

        int myid;
        int numberOfNeighs;
        int numberOfVerts;

        static utils::Logger domain_logs;

//private:
};

#endif	/* DOMAIN_H */