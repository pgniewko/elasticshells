#ifndef PAIRLIST_H
#define PAIRLIST_H

#include <vector>
#include <iostream>
#include "Vector3D.h"
#include "Box.h"
#include "Cell.h"

using namespace std;

typedef vector<int> IntList;
typedef vector<IntList> BoxList;

class Pair
{
    public:
        int k, l;
        Pair(int _k, int _l) : k(_k), l(_l) {}
};

class Pairlist
{
    public:
//    Pairlist();
        Pairlist(double mindist, int steps = 1);
//        virtual ~Pairlist();

        void compute( const  Box& dist, const Vector3D* pos, int  np, int* aindex = 0);
        static void compute_pairs(Pairlist& pl, Box& domain, vector<Cell>& cells, double r_cut, vector<Pair>& pairs);
        int findbox(Vector3D p, const Box& dist ) const;
        int num_boxes() const;
        BoxList boxlist;
        BoxList nbrlist;
        BoxList fullnbrlist;

    private:
        int n;
        double mindist;
        Vector3D margin;
        Vector3D lmin;
        Vector3D pairdist;
        int steps;
        int xb, yb, zb, xytotb;
        void set_neighbors(const Box& dist);
};
#endif