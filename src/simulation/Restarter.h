#ifndef RESTARTER_H
#define	RESTARTER_H

#include <vector>
#include <string>
#include <regex>
#include <utility>
#include <typeinfo>
#include <unordered_map>

#include "Cell.h"
#include "utils/utils.h"
#include "utils/io/XyzTraj.h"

class Restarter
{
    public:
        Restarter(std::string, std::string);
        Restarter(const Restarter& orig);
        virtual ~Restarter();

        void saveTopologyFile(const std::vector<Cell>&, std::string) const;
        void saveLastFrame(const std::vector<Cell>&) const;
        void readTopologyFile(std::vector<Cell>&) const;
        void readLastFrame(std::vector<Cell>&) const;
        void registerVMap();

    private:
        int getTotalVertices(const std::vector<Cell>&) const;

        std::pair<int, std::string> getNumberOfCells() const;
        void initCell(std::vector<Cell>&, int) const;
        void addVertices(std::vector<Cell>&, int) const;
        void addVTriangles(std::vector<Cell>&, int) const;
        void addBHinges(std::vector<Cell>&, int) const;

        //void validateVMap();

        std::string topologyFile;
        std::string lastFrameFile;
        std::unordered_map< std::string, std::pair<int, int> >  vmap;
        static utils::Logger restarter_logs;


};

#endif	/* RESTARTER_H */
