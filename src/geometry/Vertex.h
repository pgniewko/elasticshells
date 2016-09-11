#ifndef VERTEX_H
#define	VERTEX_H

#include <vector>
#include <algorithm>
#include <iostream>

#include "Environment.h"
#include "Vector3D.h"
#include "exceptions/MaxSizeException.h"
#include "exceptions/RunTimeError.h"

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
        int setCellId(int);
        int getCellId() const;

        void printVertex() const;

        void addNeighbor(int, double);
        void removeNeighbor(int);
        bool isNeighbor(int) const;
        void addTriangle(int);
        void removeTriangle(int);

        void voidForce();
        int getNumNeighs() const;
        int getNumTris() const;
        int getNeighborId(int) const;
        int getTriangleId(int) const;
        double getNeighborR0(int) const;

        Vector3D r_c;
        Vector3D f_c;

        Vector3D r_p;           // make it private
        Vector3D f_p;

        Vector3D v_p;
        Vector3D v_c;
        Vector3D a_p;
        Vector3D a_c;

        int bondedVerts[NEIGH_MAX];
        double r0[NEIGH_MAX];
        double k0[NEIGH_MAX];

        int bondedTris[TRIAN_MAX];

        int numBonded;              // make it private
        int numTris;                // make it private

        Vertex* next;


        friend std::ostream& operator<< (std::ostream&, const Vertex&);

    private:
        int myid;
        int myCellId;
};

#endif	/* VERTEX_H */