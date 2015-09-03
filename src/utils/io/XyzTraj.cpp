#include "XyzTraj.h"

XyzTraj::XyzTraj(char* tf)
{
    trajfile = tf;
}

XyzTraj::XyzTraj(const XyzTraj& orig) : trajfile(orig.trajfile) {}

XyzTraj::~XyzTraj() {}

void XyzTraj::open()
{
    os = fopen(trajfile, "w");
}

void XyzTraj::close()
{
    fclose(os);
}

//void XyzTraj::save(std::vector<Cell>& cells, int totV)
//{
//    save(cells, totV, 1.0, 1.0, 1.0);
//}

void XyzTraj::save(std::vector<Cell>& cells, int totV, double sx, double sy, double sz)
{
    int lastCellIndex = 0;
    fprintf(os, "%i\n", totV);

    for (unsigned int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            std::string strindex = new_base_index ( lastCellIndex +  cells[i].vertices[j].getId() );
            fprintf(os, "%s %10.5f %10.5f %10.5f \n", strindex.c_str(), sx * cells[i].vertices[j].xyz.x, sy * cells[i].vertices[j].xyz.y, sz * cells[i].vertices[j].xyz.z);
            
        }

        lastCellIndex += cells[i].getNumberVertices();
    }

    fflush(os);
}