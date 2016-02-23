#include "XyzTraj.h"

utils::Logger XyzTraj::xyztraj_logs("xyztraj");

XyzTraj::XyzTraj(std::string tf, std::string bf) : trajfile(tf), boxfile(bf) {}

XyzTraj::XyzTraj(const XyzTraj& orig) : trajfile(orig.trajfile), boxfile(orig.boxfile) {}

XyzTraj::~XyzTraj() {}

void XyzTraj::open()
{
    os = fopen(trajfile.c_str(), "w");

    //if ( os == NULL )
    if ( os == NULL )
    {
        xyztraj_logs << utils::LogLevel::WARNING << "Can not open file:" <<  trajfile << "\n";
    }
    
    osb = fopen(boxfile.c_str(), "w");

    if ( osb == NULL )
    {
        xyztraj_logs << utils::LogLevel::WARNING << "Can not open file:" <<  boxfile << "\n";
    }
    
}

void XyzTraj::close()
{
    if ( os != NULL )
    {
        fclose(os);
    }
    
    if ( osb != NULL )
    {
        fclose(osb);
    }
}

void XyzTraj::save(std::vector<Cell>& cells, int totV, double sx, double sy, double sz)
{
    int lastCellIndex = 0;
    fprintf(os, "%i\n", totV);

    for (unsigned int i = 0; i < cells.size(); i++)
    {
        for (int j = 0; j < cells[i].getNumberVertices(); j++)
        {
            std::string strindex = new_base_index ( lastCellIndex +  cells[i].vertices[j].getId() );
            fprintf(os, "%s %10.5f %10.5f %10.5f \n", strindex.c_str(), sx * cells[i].vertices[j].r_c.x, sy * cells[i].vertices[j].r_c.y, sz * cells[i].vertices[j].r_c.z);

        }

        lastCellIndex += cells[i].getNumberVertices();
    }

    fflush(os);
}

void XyzTraj::save_box(Box& box, double t)
{   
    fprintf(osb, "%12.6f %6.4f %6.4f %6.4f \n", t, box.getX(), box.getY(), box.getZ());   
    fflush(osb);
}