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
        XyzTraj(std::string, std::string);
        XyzTraj(const XyzTraj& orig);
        virtual ~XyzTraj();

        void open();
        void close();
        void save(const std::vector<Cell>&, int, double = 1.0, double = 1.0, double = 1.0);
        void save_box(const Box&, double);

    private:
        std::string trajfile;
        std::string boxfile;
        FILE* os;
        FILE* osb;

        static utils::Logger xyztraj_logs;
};

#endif	/* XYZTRAJ_H */

