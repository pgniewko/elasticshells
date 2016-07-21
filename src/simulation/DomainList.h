#ifndef DOMAINLIST_H
#define	DOMAINLIST_H

#include <cmath>

#include "Environment.h"
#include "geometry/Vertex.h"
#include "simulation/Box.h"
#include "Cell.h"

struct domain_t
{
    domain_t()
    {
        myid = -1;
        neighborDomainNumber = 0;
        for (int i = 0; i < MAX_D_NEIGH; i++)
            neighborDomainIdx[i] = -1;
    }
    
    int myid;
    int neighborDomainNumber;
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
        void assignVertex(Vertex*);
        int getDomainIndex(Vertex*);
        void setBoxDim(Box&);
        void setM(Box&);

        int numberofAssignedParticles();

        int getNumberOfNeigh(int);
        int getDomainNeighbor(int, int);

        double getMaxScale();
        int getVertexIdx(int, int);
        int getCellIdx(int, int);
        int getNumOfParticles(int);
        
        void calcNbForces(std::vector<Cell>&, const Box&);
        double calcNbEnergy(const std::vector<Cell>&, const Box&) const;

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

        Vertex* head[MAX_M* MAX_M* MAX_M];

        std::vector<domain_t> domains;

        static utils::Logger domainlist_logs;

        void addNeighDomain(int, int);
        void addVertex(int, int);
        
private:
    void nbForce(Vertex*, Vertex*, std::vector<Cell>&, const Box&);
    double nbEnergy(const Vertex*, const Vertex*, const std::vector<Cell>&, const Box&) const;
    bool validateLinkedDomains();

};

#endif	/* DOMAINLIST_H */