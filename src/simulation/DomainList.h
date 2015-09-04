#ifndef DOMAINLIST_H
#define	DOMAINLIST_H

#include <cmath>

#include "Environment.h"
#include "geometry/Vertex.h"
#include "simulation/Box.h"

struct domain_t
{
    int myid;
    int numberOfNeighs;
    int vertIds[MAX_IN_DOMAIN];
    int cellsIds[MAX_IN_DOMAIN];
    int neighborDomainIdx[MAX_D_NEIGH];
};

class DomainList
{
    public:
        DomainList();
        DomainList(const DomainList& orig);
        virtual ~DomainList();

        void initDomains();
        void setupDomainsList(double, Box&);
        int getDomainIndex(int, int, int);
        void voidDomains();
        void assignVertex(Vertex&, int);
        int getDomainIndex(Vertex&);
        void setBoxDim(Box&);
        void setM(Box&);

        int numberofAssignedParticles();

        int getNumberOfNeigh(int);
        int getDomainNeighbor(int, int);

        double getMaxScale();
        int getVertexIdx(int, int);
        int getCellIdx(int, int);
        int getNumOfParticles(int);

        //private:
        int m;
        int N;
        bool pbc;
        bool m_assigned;
        bool init_domains;
        double x_min, y_min, z_min;
        double x_max, y_max, z_max;
        double dx, dy, dz;
        double rc_max;

        short vertsInDomains[MAX_M* MAX_M* MAX_M];

        std::vector<domain_t> domains;

        static utils::Logger domainlist_logs;

        void addNeighDomain(int, int);
        void addVertex(int, int);

};

#endif	/* DOMAINLIST_H */