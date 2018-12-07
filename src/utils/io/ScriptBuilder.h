#ifndef SCRIPTBUILDER_H
#define	SCRIPTBUILDER_H

#include <vector>
#include <string>

#include "Shell.h"
#include "utils/utils.h"
#include "simulation/Box.h"


class ScriptBuilder
{
    public:
        ScriptBuilder(std::string, std::string, std::string, std::string);
        ScriptBuilder(const ScriptBuilder& orig);
        virtual ~ScriptBuilder();

        void save_render_script(const std::vector<Shell>&, const Box&, bool = false, double = 0.1) const;

    private:
        void print_box(const Box&, std::ofstream&) const;
        std::string script;
        std::string surfaceScript;
        std::string trajfile;
        std::string stress_script;

        static utils::Logger scriptbuilder_logs;
};

#endif	/* SCRIPTBUILDER_H */

