#include "XyzTraj.h"

XyzTraj::XyzTraj(char* tf) : names({'A','B','C','D','E','F','G','H','I','J'})
{
    trajfile = tf;
    os = fopen(trajfile, "w");
}

XyzTraj::XyzTraj(const XyzTraj& orig) {}

XyzTraj::~XyzTraj() {}

void XyzTraj::save(vector<Cell>& cells, int totV)
{
    int res1A;
    int res1B;
    int lastCellIndex = 0;
    int index;
    
    fprintf(os, "%i\n", totV);
    for (int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].numberofVertices(); j++)
        {
            index = cells[i].vertices[j].getId() + 1 + lastCellIndex;
            res1A = (int) index / 1000;
            res1B = index % 1000;
            fprintf(os, "%c%i %10.5f %10.5f %10.5f \n", names[res1A], res1B, cells[i].vertices[j].xyz.x, cells[i].vertices[j].xyz.y, cells[i].vertices[j].xyz.z);
        }
        
        lastCellIndex += cells[i].numberofVertices();
    }
    
}
void XyzTraj::close()
{
    fclose(os);
}

//void Simulator::write_pos_traj(const char* filename, bool wrap)
//{
//     ofstream os(filename, ios::app);


