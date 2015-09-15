#ifndef SCRIPTBUILDER_H
#define	SCRIPTBUILDER_H

#include <vector>

#include "Cell.h"
#include "utils/utils.h"
#include "simulation/Box.h"


class ScriptBuilder
{
    public:
        ScriptBuilder(char*, char*, char*, char*);
        ScriptBuilder(const ScriptBuilder& orig);
        virtual ~ScriptBuilder();

        void printBox(std::ofstream&, const Box&);
        void saveRenderScript(const std::vector<Cell>&, const Box&, bool = false, double = 0.1);
        void saveSurfaceScript(const std::vector<Cell>&);

        void saveStressScript(std::vector<Cell>&, const Box&);
        void saveBfactors();
        void setDrawBox(bool);
    private:
        char* script;
        char* stress_script;
        char* surfaceScript;
        char* trajfile;
        bool drawBox;
};

#endif	/* SCRIPTBUILDER_H */

