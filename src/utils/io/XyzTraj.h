#ifndef XYZTRAJ_H
#define	XYZTRAJ_H

#include <fstream>
#include <vector>

#include "Cell.h"
#include "Environment.h"
#include "simulation/Box.h"
#include "../utils.h"

class XyzTraj
{
    public:
        XyzTraj(char*);
        XyzTraj(const XyzTraj& orig);
        virtual ~XyzTraj();

        void open();
        void close();
//        void save(std::vector<Cell>&, int);
        void save(std::vector<Cell>&, int, double=1.0, double=1.0, double=1.0);

    private:
        
        char* trajfile;
        FILE* os;
};

#endif	/* XYZTRAJ_H */

