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
    
    bool operator > (const nblist_t& rhs) 
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

    bool operator < (const nblist_t& rhs) 
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

class Vertex
{
    public:
        Vertex();
        Vertex(double, double, double);
        Vertex(Vector3D);
        Vertex(const Vertex& orig);
        virtual ~Vertex();
        int setId(int);
        int getId();
        double setMass(double);
        double getMass();
        double setVisc(double);
        double getVisc();

        void printVertex();

        void addNeighbor(int, double);
        void addNbNeighbor(int, int);
        bool isNeighbor(int);
        void addTriangle(int);

        void voidForce();
        void voidVelocity();

        int getNumNeighs();
        int getNumTris();
        int getNeighborId(int);
        int getTriangleId(int);
        double getNeighborR0(int);
        
        void sortNbList();

        Vector3D xyz;
        Vector3D force;
        Vector3D velocity;

        Vector3D tmp_xyz;           // make it private
        Vector3D tmp_force;         // make it private
        Vector3D tmp_velocity;      // make it private

        int bondedVerts[NEIGH_MAX];
        int bondedTris[TRIAN_MAX];
        double r0[NEIGH_MAX];

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
};

#endif	/* VERTEX_H */