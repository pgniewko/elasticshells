#ifndef XYZTRAJ_H
#define	XYZTRAJ_H

#include <vector>
#include <string>

#include "Shell.h"
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
        void open_traj();
        void open_lf();
        void open_box();
        void close();
        void close_traj();
        void close_box();

        void save_traj(const std::vector<Shell>&, int);
        void save_traj(const std::vector<Shell>&, int, const Box&);
        void save_box(const Box&, double);

        const std::vector<std::string> read_saved_box() const;
        uint count_frames() const;
        const std::string get_traj_file() const
        {
            return trajfile;
        }

    private:
        std::string trajfile;
        std::string boxfile;
        FILE* os;
        FILE* osb;

        static utils::Logger xyztraj_logs;
};

#endif	/* XYZTRAJ_H */

