#ifndef RESTARTER_H
#define	RESTARTER_H

#include <vector>
#include <string>
#include <utility>
#include <typeinfo>
#include <unordered_map>

#include "Shell.h"
#include "geometry/Vertex.h"
#include "utils/utils.h"
#include "utils/io/XyzTraj.h"

class Restarter
{
    public:
        Restarter(std::string, std::string);
        Restarter(const Restarter& orig);
        virtual ~Restarter();

        void saveTopologyFile(const std::vector<Shell>&, std::string) const;
        void saveLastFrame(const std::vector<Shell>&) const;
        void readTopologyFile(std::vector<Shell>&) const;
        void readLastFrame(std::vector<Shell>&) const;
        void registerVMap();

        void readFrame(std::string, std::vector<Shell>&, int) const;
        void assignTurgors(std::string, std::vector<Shell>&) const;
        void assignBoxSize(std::string, Box&) const;

    private:
        int getTotalVertices(const std::vector<Shell>&) const;

        std::pair<int, std::string> getNumberOfCells() const;
        void initCell(std::vector<Shell>&, int) const;
        void addVertices(std::vector<Shell>&, int) const;
        void addVTriangles(std::vector<Shell>&, int) const;
        void addBHinges(std::vector<Shell>&, int) const;

        //void validateVMap();

        std::string topologyFile;
        std::string lastFrameFile;
        std::unordered_map< std::string, std::pair<int, int> >  vmap;
        static utils::Logger restarter_logs;


};

#endif	/* RESTARTER_H */

