#ifndef SCRIPTBUILDER_H
#define	SCRIPTBUILDER_H

#include <fstream>
#include <vector>

#include "Cell.h"
#include "simulation/Box.h"

class ScriptBuilder
{
    public:
        ScriptBuilder(char*, char*, char*);
        ScriptBuilder(const ScriptBuilder& orig);
        virtual ~ScriptBuilder();

        void printBox(std::ofstream&, Box&);
        void saveRenderScript(std::vector<Cell>&, Box&, bool box = false);
        void saveSurfaceScript(std::vector<Cell>&);

        void setDrawBox(bool);
    private:

        //char names[10];
        char* script;
        char* surfaceScript;
        char* trajfile;
        bool drawBox;

};

#endif	/* SCRIPTBUILDER_H */

