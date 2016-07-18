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
    //int vertIds[MAX_IN_DOMAIN];
    //int cellsIds[MAX_IN_DOMAIN];
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

        //short vertsInDomains[MAX_M* MAX_M* MAX_M];
        Vertex* head[MAX_M* MAX_M* MAX_M];

        std::vector<domain_t> domains;

        static utils::Logger domainlist_logs;

        void addNeighDomain(int, int);
        void addVertex(int, int);
        
        
        inline void newGetDistance(Vector3D& dkj, Vector3D& vj, Vector3D& vk, const Box& box) const
        {
        dkj = vk - vj;

        if (box.pbc)
        {
            double x, y, z;
            double bsx = 2 * box.getX();
            double bsy = 2 * box.getY();
            double bsz = 2 * box.getZ();
            x = round(dkj.x / bsx) *  bsx;
            y = round(dkj.y / bsy) *  bsy;
            z = round(dkj.z / bsz) *  bsz;
            dkj.x -= x;
            dkj.y -= y;
            dkj.z -= z;
        }
    }
        
private:
    void nbForce(Vertex*, Vertex*, std::vector<Cell>&, const Box&);
    bool validateLinkedDomains();

};

#endif	/* DOMAINLIST_H */