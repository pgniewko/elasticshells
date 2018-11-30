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
        Vertex(const Vertex& orig);
        virtual ~Vertex();
        int set_id(int);
        int get_id() const;
        int set_shell_id(int);
        int get_shell_id() const;

        void printVertex() const;

        void addNeighbor(int);
        bool isNeighbor(int) const;
        void addTriangle(int);

        void voidForce();
        int get_vertex_degree() const;
        int getNumTris() const;
        int getNeighborId(int) const;
        int getTriangleId(int) const;

        Vector3D r_c;
        Vector3D f_c;

        Vector3D r_p;           // make it private
        Vector3D f_p;

        Vector3D v_p;
        Vector3D v_c;
        Vector3D a_p;
        Vector3D a_c;

        std::vector<int> bondedVerts;
        std::vector<int> bondedTris;

        int vertex_degree;              // make it private
        int facets_number;                // make it private

        Vertex* next = NULL;


        friend std::ostream& operator<< (std::ostream&, const Vertex&);

    private:
        int myid;
        int myCellId;
};

#endif	/* VERTEX_H */