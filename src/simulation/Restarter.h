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

        void save_topology_file(const std::vector<Shell>&) const;
        void save_last_frame(const std::vector<Shell>&, const Box&) const;
        void read_topology_file(std::vector<Shell>&) const;
        void read_last_frame(std::vector<Shell>&) const;
        void register_vmap();

        void read_frame(std::string, std::vector<Shell>&, int) const;
        void assign_turgors(std::string, std::vector<Shell>&) const;
        void assign_box_size(std::string, Box&) const;

    private:
        int get_total_vertices(const std::vector<Shell>&) const;

        std::pair<int, std::string> get_number_of_shells() const;
        void init_shell(std::vector<Shell>&, int) const;
        void add_vertices(std::vector<Shell>&, int) const;
        void add_elements(std::vector<Shell>&, int) const;
        void add_hinges(std::vector<Shell>&, int) const;

        std::string topologyFile;
        std::string lastFrameFile;
        std::unordered_map< std::string, std::pair<int, int> >  vmap;
        static utils::Logger restarter_logs;
};

#endif	/* RESTARTER_H */

