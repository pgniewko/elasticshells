#ifndef XYZTRAJ_H
#define	XYZTRAJ_H

#include <fstream>
#include <vector>

#include "Cell.h"
#include "Environment.h"
#include "simulation/Box.h"

class XyzTraj
{
    public:
        XyzTraj(char*);
        XyzTraj(const XyzTraj& orig);
        virtual ~XyzTraj();

        void open();
        void close();
        void save(std::vector<Cell>&, int);
        void save(std::vector<Cell>&, int, double, double, double);

    private:
        char* trajfile;
        FILE* os;
        //char names[10];
};

#endif	/* XYZTRAJ_H */

