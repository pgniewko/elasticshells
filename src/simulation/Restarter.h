#ifndef RESTARTER_H
#define	RESTARTER_H

#include <vector>
#include <string>

#include "Cell.h"
#include "utils/utils.h"
#include "utils/io/XyzTraj.h"

class Restarter
{
    public:
        Restarter(std::string, std::string);
        Restarter(const Restarter& orig);
        virtual ~Restarter();

        void saveTopologyFile(const std::vector<Cell>&, char*) const;
        void saveLastFrame(const std::vector<Cell>&) const;
        void readTopologyFile(const std::vector<Cell>&) const;
        void readLastFrame(const std::vector<Cell>&) const;

    private:
        int getTotalVertices(const std::vector<Cell>&) const;

        std::string topologyFile;
        std::string lastFrameFile;
        static utils::Logger restarter_logs;

};

#endif	/* RESTARTER_H */

