#ifndef XYZTRAJ_H
#define	XYZTRAJ_H


#include <fstream>
#include <float.h>      /* DBL_MAX */
#include <cstring>
#include <vector>
#include <stdio.h>      /* fprintf*/

#include "Cell.h"
#include "simulation/Box.h"

class XyzTraj {
public:
    XyzTraj(char *);
    XyzTraj(const XyzTraj& orig);
    virtual ~XyzTraj();
    
    void save(vector<Cell>&, int);
    void close();
private:
    char* trajfile;
    FILE* os;
    char names[10];
};

#endif	/* XYZTRAJ_H */

