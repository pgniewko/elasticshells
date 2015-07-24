#ifndef SCRIPTBUILDER_H
#define	SCRIPTBUILDER_H

#include <fstream>
#include <vector>

#include "Cell.h"
#include "simulation/Box.h"

class ScriptBuilder
{
    public:
        ScriptBuilder(char*, char*, char*, char*);
        ScriptBuilder(const ScriptBuilder& orig);
        virtual ~ScriptBuilder();

        void printBox(std::ofstream&, Box&);
        void saveRenderScript(std::vector<Cell>&, Box&, bool = false, double = 0.1);
        void saveSurfaceScript(std::vector<Cell>&);

        void saveStressScript(std::vector<Cell>&, Box&);
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

