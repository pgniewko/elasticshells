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

        void saveRenderScript(const std::vector<Cell>&, const Box&, bool = false, double = 0.1);
        void saveSurfaceScript(const std::vector<Cell>&);

        void saveStressScript(const std::vector<Cell>&, const Box&);
        void saveBfactors();
        void setDrawBox(bool);
    private:
        void printBox(std::ofstream&, const Box&);
        //char* script;
        //char* stress_script;
        //char* surfaceScript;
        //char* trajfile;
        
        std::string script;
        std::string surfaceScript;
        std::string trajfile;
        std::string stress_script;
        bool drawBox;
        
        static utils::Logger scriptbuilder_logs;
};

#endif	/* SCRIPTBUILDER_H */

