#ifndef XYZTRAJ_H
#define	XYZTRAJ_H

#include <vector>
#include <string>

#include "Cell.h"
#include "Environment.h"
#include "simulation/Box.h"
#include "utils/utils.h"

class XyzTraj
{
    public:
        XyzTraj(std::string);
        XyzTraj(const XyzTraj& orig);
        virtual ~XyzTraj();

        void open();
        void close();
        void save(std::vector<Cell>&, int, double = 1.0, double = 1.0, double = 1.0);

    private:
        std::string trajfile;
        FILE* os;

        static utils::Logger xyztraj_logs;
};

#endif	/* XYZTRAJ_H */

