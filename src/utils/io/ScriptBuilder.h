#ifndef SCRIPTBUILDER_H
#define	SCRIPTBUILDER_H

#include <fstream>
#include <float.h>      /* DBL_MAX */
#include <cstring>
#include <vector>
#include <stdio.h>      /* fprintf*/

#include "Cell.h"
#include "simulation/Box.h"

class ScriptBuilder {
public:
    ScriptBuilder(char*, char*, char*);
    ScriptBuilder(const ScriptBuilder& orig);
    virtual ~ScriptBuilder();
    
    void printBox(ofstream&, Box&);
    void saveRenderScript(vector<Cell>&, Box&, bool box=false);
    void saveSurfaceScript(vector<Cell>&);
    
    //void setRenderScriptName();
    //void setSurfScriptName();
    void setDrawBox(bool);
private:

    char names[10];// = {'A','B','C','D','E','F','G','H','I','J'};
    char* script;
    char* surfaceScript;
    char* trajfile;
    bool drawBox;
    
};

#endif	/* SCRIPTBUILDER_H */

