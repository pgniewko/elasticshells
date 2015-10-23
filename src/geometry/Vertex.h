#ifndef VERTEX_H
#define	VERTEX_H

#include <vector>
#include <algorithm>
#include <iostream>

#include "Environment.h"
#include "Vector3D.h"
#include "exceptions/MaxSizeException.h"
#include "exceptions/RunTimeError.h"

struct nblist_t
{
    int cell_id;
    int vertex_id;

    bool operator > (const nblist_t& rhs) const
    {
        if (cell_id == rhs.cell_id)
        {
            return vertex_id > rhs.vertex_id;
        }
        else
        {
            return cell_id > rhs.cell_id;
        }
    }

    bool operator < (const nblist_t& rhs) const
    {
        if (cell_id == rhs.cell_id)
        {
            return vertex_id < rhs.vertex_id;
        }
        else
        {
            return cell_id < rhs.cell_id;
        }
    }
};

enum class vertex_t
{
    MOTHER,
    BUD,
    GHOST
};

class Vertex
{
        friend class Tinker;
    public:
        Vertex();
        Vertex(double, double, double);
        Vertex(Vector3D);
        Vertex(const Vertex& orig);
        virtual ~Vertex();
        int setId(int);
        int getId() const;
        double setMass(double);
        double getMass() const;
        double setVisc(double);
        double getVisc() const;

        void printVertex();

        void addNeighbor(int, double);
        void removeNeighbor(int);
        void addNbNeighbor(int, int);
        bool isNeighbor(int) const;
        void addTriangle(int);
        void removeTriangle(int);

        void voidForce();
        void voidVelocity();

        int getNumNeighs() const;
        int getNumTris() const;
        int getNeighborId(int) const;
        int getTriangleId(int) const;
        double getNeighborR0(int) const;

        void sortNbList();
        void normalizedR0(double);

        void addTime(double);
        void voidTime();

        const vertex_t& getMyType() const;

        Vector3D r_c;
        Vector3D v_c;
        Vector3D f_c;
        
        Vector3D r_p;           // make it private
        Vector3D v_p;           // make it private
        Vector3D f_p;
        
        
        int bondedVerts[NEIGH_MAX];
        double r0[NEIGH_MAX];

        int bondedTris[TRIAN_MAX];


        int nbVerts[NBNEI_MAX];
        int nbCellsIdx[NBNEI_MAX];

        int numBonded;
        int numTris;
        int numNbNeighs;

        int domainIdx;              // make it private

    private:
        int myid;
        double mass;
        double visc;

        double gtimer;
        vertex_t my_type;
};

#endif	/* VERTEX_H */