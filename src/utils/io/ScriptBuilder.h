#ifndef SCRIPTBUILDER_H
#define	SCRIPTBUILDER_H

#include <vector>
#include <string>

#include "Cell.h"
#include "utils/utils.h"
#include "simulation/Box.h"


class ScriptBuilder
{
    public:
        ScriptBuilder(std::string, std::string, std::string, std::string);
        ScriptBuilder(const ScriptBuilder& orig);
        virtual ~ScriptBuilder();

        void saveRenderScript(const std::vector<Cell>&, const Box&, bool = false, double = 0.1) const;
        void saveSurfaceScript(const std::vector<Cell>&) const;
        void saveContactStress(const std::vector<Cell>&, const Box&) const;
        void saveStrainScript(const std::vector<Cell>&, const Box&) const;

    private:
        void printBox(const Box&, std::ofstream&) const;
        std::string script;
        std::string surfaceScript;
        std::string trajfile;
        std::string stress_script;

        static utils::Logger scriptbuilder_logs;
};

#endif	/* SCRIPTBUILDER_H */

